
fin0, sdir, fout = ARGV
fout = open(fout, "w")
# [1-6: BED6] OM-RGC.v1.000702001     409     695     M02178:13:000000000-A9RUR:1:2107:19743:8499/2   2       -       
# [7-12: fin] low_score scr:256/457.60 len:286/(150.5&50) strand:OK distant:OK 2S286M13S
# [13- : add] type:[ABCX] weight:[1, 0.5, 0] min_aln_len:([150.5, 50].max) read_len:(2*150.5) 

type2weight   = %w|A B C X|.zip([0.5, 1,   1, 0]).inject({}){ |h, (i, j)| h[i] = j; h }
type2num_read = %w|A B C X|.zip([0.5, 0.5, 1, 0]).inject({}){ |h, (i, j)| h[i] = j; h }
# A: pe,       same      => weight:0.5
# B: pe,       different => weight:1
# C: pe_or_up, only_one  => weight:1
# X: others              => weight:0

# calc by total num read
# not  by sum_weight

# parse reference seq length
ref2len = {}
IO.read(fin0).split(/^>/)[1..-1].each{ |ent|
	lab, *seq = ent.split("\n")
	lab = lab.split[0]
	len = seq.join.gsub(/\s+/, "").size
	ref2len[lab] = len
}

# data stores
ref2types               = Hash.new{ |h, i| h[i] = Hash.new(0) }
ref2ref_sum_weight      = Hash.new(0)
ref2ref_weight_per_kb   = Hash.new(0) # sum of (weight * 1000 / (ref_len + read_len - 2 * min_aln_len)) of each mapping
total_num_read          = 0
total_sum_weight        = 0

# parse mapping info
Dir["#{sdir}/*.bed.classified"].each{ |fin|
	IO.readlines(fin).each{ |l|
		ref, type, weight, min_aln_len, read_len = l.chomp.split("\t").values_at(0, 10..13)

		# count type
		type        = type.split(":")[1]
		ref2types[ref][type] += 1
		next if type == "X"

		# parse weight, length, ,,,
		ref_len = ref2len[ref]
		raise unless ref_len

		# weight      = (weight.split(":")[1]      || false).to_f
		weight      = type2weight[type].to_f
		num_read    = type2num_read[type]
		read_len    = (read_len.split(":")[1]    || false).to_f
		min_aln_len = (min_aln_len.split(":")[1] || false).to_f

		# add weight (for each ref, for all)
		ref2ref_sum_weight[ref] += weight
		total_sum_weight        += weight
		total_num_read          += num_read

		# calc weighted_count
		# [TODO] make option for advanced num_frame mode
		# num_frame       = ref_len + read_len - 2 * min_aln_len
		# (p [ref_len, read_len, min_aln_len, num_frame]; puts l; raise) if num_frame < 0
		num_frame     = ref_len
		weight_per_kb = weight * 1000 / num_frame

		# p [weight_per_kb, weight, num_frame, ref_len, read_len, 2 * min_aln_len]

		ref2ref_weight_per_kb[ref] += weight_per_kb

	}
}
p [sdir, total_num_read, total_sum_weight]

# def calc_FPKM(weight_per_kb, total_sum_weight)
# 	return weight_per_kb * 1_000_000 / total_sum_weight
# end
def calc_FPKM(weight_per_kb, total_num_read)
	return weight_per_kb * 1_000_000 / total_num_read
end

header = %w|ID length FPKM #mapped mapping_type|
fout.puts header*"\t"

ref2len.each{ |ref, ref_len|
	# types
	types = ref2types[ref]
	type_str = (types.size != 0) ? types.keys.sort.map{ |type| [type, types[type]]*":" }*"," : "-"

	# fpkm
	weight_per_kb = ref2ref_weight_per_kb[ref]

	# fpkm = calc_FPKM(weight_per_kb, total_sum_weight)
	fpkm = calc_FPKM(weight_per_kb, total_num_read)
	fpkm_str = "%.6f" % fpkm

	# sum_weight_for_ref
	sum_weight = ref2ref_sum_weight[ref].to_i
	# sum_weight = ref2ref_sum_weight[ref]
	# sum_weight_str = "%.6f" % sum_weight

	# fout.puts [ref, ref_len, fpkm_str, sum_weight_str, type_str]*"\t"
	fout.puts [ref, ref_len, fpkm_str, sum_weight, type_str]*"\t"
}
