#!/bin/bash
#
#  CountMappedReads.rake - a command line tool for mapping/counting NGS reads on reference genomes or genes
#
#    Copyright: 2017 (C) Yosuke Nishimura (yosuke@kuicr.kyoto-u.ac.jp)
#    License: MIT license
#    Initial version: 2017-03-03
#


# {{{ procedures
WriteBatch  = lambda do |outs, jdir, t|
	outs.each_slice(10000).with_index(1){ |ls, idx|
		open("#{jdir}/#{t.name.split(".")[1..-1]*"."}.#{idx}", "w"){ |fjob|
			fjob.puts ls
		}
	}
end

RunBatch    = lambda do |jdir, queue, nthreads, mem, wtime, ncpus|
	# [TODO] queue validation
	Dir["#{jdir}/*"].sort_by{ |fin| fin.split(".")[-1].to_i }.each{ |fin|
		if queue != ""
			raise("`--queue #{queue}': invalid queue") if %w|JP1 cdb|.include?(queue)
			sh "qsubarraywww -q #{queue} -l ncpus=#{nthreads} -l mem=#{mem}gb -l walltime=#{wtime} #{fin}"
		elsif ncpus != ""
			raise("`--ncpus #{ncpus}': not an integer") if ncpus !~ /^\d+$/
			sh "parallel --jobs #{ncpus} <#{fin}"
		else
			sh "sh #{fin}"
		end
	}
end

PrintStatus = lambda do |current, total, status, t|
	puts ""
	puts "===== #{Time.now}"
	puts "===== step #{current} / #{total} (#{t.name}) -- #{status}"
	puts ""
	$stdout.flush
end
# }}} procedures


# {{{ default (run all tasks)
task :default do
	tasks = %w|
	01-1.bowtie2_build 01-2.bowtie2 01-3.samtools
	02-1.mapping_selection 02-2.classify_mapping_type_A-C 02-3.make_list_of_FPKM
	|
	NumStep = tasks.size
	tasks.each.with_index(1){ |task, idx|
		# next if task =~ /^03-2/ and ENV["notree"] != "" # skip tree calculation if notree mode
		Rake::Task[task].invoke(idx)
	}
	# begin
	# rescue SystemExit => e
	# 	puts e.status
	# end
end
# }}} default (run all tasks)


# {{{ tasks 01
desc "01-1.bowtie2_build"
task "01-1.bowtie2_build", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	dir      = ENV["dir"]
	fin      = ENV["fin"]
	threads  = ENV["threads"]||"1"
	idir     = "#{dir}/bowtie2-index"
	fa       = "#{idir}/#{File.basename(fin)}"

	mkdir_p idir
	sh "cp #{fin} #{fa}"
	sh "bowtie2-build --threads #{threads} #{fa} #{fa}"
end
desc "01-2.bowtie2"
task "01-2.bowtie2", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	dir       = ENV["dir"]
	fin       = ENV["fin"]
	threads   = ENV["threads"]||"1"
	pe1       = ENV["pe1"]||""
	pe2       = ENV["pe2"]||""
	up        = ENV["up"]||""
	mapping   = ENV["mapping"]||"end-to-end" # "end-to-end" or "local"
	score_min = ENV["score_min"]||""
	score_min = (mapping == "end-to-end" ? "L,0,-0.6" : "L,0,0.6") if score_min == ""
	# end-to-end, ~90% idt: "L,0,-0.6"
	# end-to-end, ~95% idt: "L,0,-0.3"
	# local, ~90% idt: "L,0,0.6" # initial screening (for 50% length alignment)
	# local, ~95% idt: "L,0,0.8" # initial screening (for 50% length alignment)
	idir      = "#{dir}/bowtie2-index"
	fa        = "#{idir}/#{File.basename(fin)}"
	mdir      = "#{dir}/mapping"
	sam       = "#{mdir}/out.sam"

	input     = ""
	input    += "-1 #{pe1} " if pe1 != ""
	input    += "-2 #{pe2} " if pe2 != ""
	input    += "-U #{up} "  if up  != ""

	mkdir_p mdir
	sh "bowtie2 --#{mapping} --score-min #{score_min} --threads #{threads} -x #{fa} #{input} -S #{sam}"
end
desc "01-3.samtools"
task "01-3.samtools", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	dir       = ENV["dir"]
	threads   = ENV["threads"]||"1"
	mem       = ENV["mem"]||"1G"
	mdir      = "#{dir}/mapping"
	sam       = "#{mdir}/out.sam"

	bam       = sam.gsub(/\.sam$/, ".bam.tmp")
	sorted    = sam.gsub(/\.sam$/, ".bam")
	map_q0    = sam.gsub(/\.sam$/, ".mapped") # mapped with Q > 0

	str =  ["samtools view --threads #{threads} -bS #{sam} -o #{bam}",
				  "rm #{sam}",
	        "samtools sort --threads #{threads} -m #{mem} -O bam -T #{mdir} #{bam} -o #{sorted}",
				  "rm #{bam}",
	        "samtools index #{sorted}",
	        "samtools view -F 0x4 #{sorted} >#{map_q0}"].join(" && ")

	sh "#{str}"
end
# }}} tasks 01


# {{{ tasks 02
desc "02-1.mapping_selection"
task "02-1.mapping_selection", ["step"] do |t, args|
	dir       = ENV["dir"]
	mapping   = ENV["mapping"]||"end-to-end" # "end-to-end" or "local"
	mdir      = "#{dir}/mapping"
	map_q0    = "#{mdir}/out.mapped"
	sdir      = "#{dir}/mapping/selection"

	script    = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"

	mkdir_p sdir
	sh "ruby #{script} #{map_q0} #{sdir} #{mapping}"
end
desc "02-2.classify_mapping_type_A-C"
task "02-2.classify_mapping_type_A-C", ["step"] do |t, args|
	dir       = ENV["dir"]
	sdir      = "#{dir}/mapping/selection"

	script    = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"

	Dir["#{sdir}/*.bed"].each{ |bed|
		fout = "#{bed}.classified"
		sh "ruby #{script} #{bed} #{fout}"
	}
end
desc "02-3.make_list_of_FPKM"
task "02-3.make_list_of_FPKM", ["step"] do |t, args|
	dir       = ENV["dir"]
	fin       = ENV["fin"]
	sdir      = "#{dir}/mapping/selection"

	rdir      = "#{dir}/result"
	ffpkm     = "#{rdir}/count-FPKM.tsv"

	mkdir_p rdir
	script    = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"

	sh "ruby #{script} #{fin} #{sdir} #{ffpkm}"
end
# }}} tasks 02

