
fin, odir, mode = ARGV
#   fin: SAM of mapping (>=q0)
# fstat: stats
# fout:  lines of each type
# mode: end-to-end or local

stat = Hash.new{ |h, i| h[i] = {} }
fout = Hash.new{ |h, i| h[i] = {} }
#stat[<pe_or_me>][<mapping>][<state>] == count
stat["pe"]["same"]      = %w|ok too_short low_score both|.inject(Hash.new(0)){ |h, i| h[i] = 0; h }
stat["pe"]["different"] = %w|ok too_short low_score both|.inject(Hash.new(0)){ |h, i| h[i] = 0; h }
stat["pe"]["only_one"]  = %w|ok too_short low_score both|.inject(Hash.new(0)){ |h, i| h[i] = 0; h }
stat["up"]["only_one"]  = %w|ok too_short low_score both|.inject(Hash.new(0)){ |h, i| h[i] = 0; h }

stat.each{ |pe_or_up, h1|
	h1.each{ |mapping, h2|
		fout[pe_or_up][mapping]  = %w|sam bed|.inject({}){ |h, i| h[i] = open("#{odir}/#{pe_or_up}-#{mapping}.#{i}", "w"); h }
	}
}

fstat = open("#{odir}/mapping.stat", "w")

IO.readlines(fin).each{ |l|
	# flags
	both_mapped    = false #     for paired end mapping
	wrong_stranded = false #     for judge concordant / discordant
	too_short_aln  = false # [1] filter out unless (aln_read >= 50 and aln_read >= read_len * 0.5)
	too_low_score  = false # [2] filter out unless (aln_score >= 1.2 * aln_read)
	too_distant    = false #     for judge concordant / discordant
	pe_or_up       = false #     pe: paired-end read, up: unpaired read
	first_or_last  = false #     first: in 1st file, last: in 2nd file
	compl_mapping  = false #     read is complementally mapped ?

	## parse SAM line
	qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, as_i = l.chomp.split(/\t/).values_at(0..8, 11)
	_cigar = cigar.dup

	## parse flag
	flag = "%012d" % flag.to_i.to_s(2)
	# (1) 0x1: template having multiple segments in sequencing => 1: paired, 0: single
	# (7) 0x40 the first segment in the template
	# (8) 0x80 the last segment in the template
	pe_or_up = flag[-1] == "1" ? "pe" : "up"
	if pe_or_up == "pe"
		first_or_last = flag[-7] == "1" ? "/1" : (flag[-8] == "1" ? "/2" : raise)
	end

	# (3) 0x4 segment unmapped
	# (4) 0x8 next segment in the template unmapped
	both_mapped = true if pe_or_up == "pe" and flag[-3] == "0" and flag[-4] == "0"
	# (5) 0x10 SEQ being reverse complemented
	# (6) 0x20 SEQ of the next segment in the template being reversed
	# wrong_stranded = true if both_mapped and flag[-5] == flag[-6] # [!!!] disabled
	compl_mapping  = true if flag[-5] == "1"

	## parse cigar
	# {"S"=>7, "M"=>291, "I"=>1, "D"=>1}
	cigar = cigar.split(/([A-Z])/).each_slice(2).inject(Hash.new(0)){ |h, i| h[i[1]] += i[0].to_i; h }

	# check inclusion of the other cigar symbol (N, H, P, =, X)
	raise if cigar.keys.join('') =~ /[^MDIS]/ 

	# calc length
	read_len = %w|M I S|.inject(0){ |s, i| s += cigar[i] }
	refr_len = %w|M D  |.inject(0){ |s, i| s += cigar[i] }
	aln_real = %w|M I D|.inject(0){ |s, i| s += cigar[i] }
	aln_read = %w|M I  |.inject(0){ |s, i| s += cigar[i] }

	# [1] check too short alignment if local mode
	if mode == "local" and (aln_read < 50 or aln_read < read_len * 0.5)
		too_short_aln = true
	end

	## parse mapping score
	aln_score = as_i[/^AS:i:(\d+)$/, 1].to_i

	# [2] check too low score of alignment if local mode
	# min_score = 1.6 * aln_real # 95% identity if gap = 0
	# min_score = 1.6 * aln_read # 95% identity if gap = 0
	min_score = 1.2 * aln_read # 90% identity if gap = 0
	if mode == "local" and aln_score < min_score
		too_low_score = true
	end


	## parse tlen and validate ( -1000 <= tlen <= 1000 )
	# tlen = tlen.to_i
	# too_distant = true if tlen < -1000 or tlen > 1000 # [!!!] disabled

	## parse rnext
	mapping = if both_mapped
							case rnext
							when "=" # paired concordant mapping
								# (wrong_stranded or too_distant) ? "discordant" : "concordant"
								"same"
							when "*" then (puts l;raise)
							else # paired mapping on different sequences
								"different"
							end
						else "only_one" # single read or only one side mapping of PE
						end

	## state
	state = if too_short_aln and too_low_score then "both"
					elsif too_short_aln then "too_short"
					elsif too_low_score then "low_score"
					else "ok"
					end

	## make output of stats
	stat[pe_or_up][mapping][state] += 1

	## make bed output
	query       = pe_or_up == "up" ? qname : (qname + first_or_last)
	start       = pos.to_i - 1 # 1-based => 0-based
	stop        = start + refr_len
	strand      = compl_mapping ? "-" : "+"
	score_str   = "scr:%d/%.2f"    % [aln_score, min_score]
	alnlen_str  = "len:%d/(%.1f&%d)" % [aln_read, read_len * 0.5, 50]
	# strand_str  = (wrong_stranded and rnext == "=") ? "strand:NG"  : "strand:OK"
	# distant_str = (too_distant    and rnext == "=") ? "distant:NG" : "distant:OK"
	# bed = [rname, start, stop, query, mapq, strand, state, score_str, alnlen_str, strand_str, distant_str, _cigar]
	bed = [rname, start, stop, query, mapq, strand, state, score_str, alnlen_str, _cigar]

	## make output
	fout[pe_or_up][mapping]["sam"].puts l
	fout[pe_or_up][mapping]["bed"].puts bed*"\t"
}

header = %w|pe_or_up mapping OK TooShort LowScore TooShort&LowScore|
fstat.puts header*"\t"
stat.each{ |pe_or_up, h|
	h.each{ |mapping, g|
		counts = %w|ok too_short low_score both|.map{ |state| g[state] }
		fstat.puts [pe_or_up, mapping, counts]*"\t"
	}
}
