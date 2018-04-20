
fin, fout = ARGV

fout = open(fout, "w")
# [1-6: BED6] OM-RGC.v1.000702001     409     695     M02178:13:000000000-A9RUR:1:2107:19743:8499/2   2       -       
# [7-10: fin] low_score scr:256/457.60 len:286/(150.5&50) 2S286M13S
# [11- : add] type:[ABCX] weight:[1, 0.5, 0] min_aln_len:([150.5, 50].max) read_len:(2*150.5) 

type2weight = %w|A B C X|.zip([0.5, 1, 1, 0]).inject({}){ |h, (i, j)| h[i] = j; h }
# A: pe,       same      => weight:0.5
# B: pe,       different => weight:1
# C: pe_or_up, only_one  => weight:1
# X: others              => weight:0

## mapping state
# selfstate: ok or not
# nextstate: ok or not (the other of pair) for concordant and discordant

# fetch nextstate
case fin
when /pe-same\.bed$/, /pe-different\.bed$/ # pe-same or pe-different
	readpair2linepair = Hash.new{ |h, i| h[i] = [] }
	read2nextstate    = {} # h["read/1"] => <state of read/2>, h["read/2"] => <state of read/1>
	IO.readlines(fin).each{ |l|
		readpair = l.chomp.split(/\t/)[3][/^([^\/]+)\/[12]/, 1]
		raise unless readpair
		readpair2linepair[readpair] << l
	}
	readpair2linepair.each{ |readpair, lines|
		raise if lines.size != 2
		reads  = lines.map{ |l| l.split(/\t/)[3] }
		states = lines.map{ |l| l.split(/\t/)[6] }.reverse
		reads.zip(states){ |read, state| read2nextstate[read] = state }
	}
end

# parse fin and add additional columns
IO.readlines(fin).each{ |l|
	a = l.chomp.split(/\t/)

	read, selfstate = a.values_at(3, 6)
	read_half_len   = (a[8][/len:\d+\/\((\d+\.\d)&50\)/, 1] || false).to_f

	min_aln_len_str = "min_aln_len:%.1f" % [read_half_len, 50].max
	read_len_str    = "read_len:%d"      % (2 * read_half_len)

	# assign type
	type = case fin
				 when /-only_one\.bed$/ #only-one
					 selfstate == "ok" ? "C" : "X"
				 when /pe-different\.bed$/ # pe-different
					 nextstate = read2nextstate[read]
					 raise unless nextstate
					 selfstate == "ok" ? (nextstate == "ok" ? "B" : "C") : "X"
				 when /pe-same\.bed$/ # pe-same
					 nextstate = read2nextstate[read]
					 raise unless nextstate
					 selfstate == "ok" ? (nextstate == "ok" ? "A" : "C") : "X"
				 end

	# output
	type_str   = "type:#{type}"
	weight_str = "weight:%.1f" % type2weight[type]
	fout.puts [a, type_str, weight_str, min_aln_len_str, read_len_str]*"\t"
}
