import sys
import os

genome_file="E_faecalis.fasta"
output_file="test_TA.tab"
patt = "TA"


header_done=False
pos = 0
last_nuc=""

out_handle = open(output_file, 'w')
with open(genome_file) as in_handle:
	while True:
		c = in_handle.read(1)
		if not c:
			print("End of file")
			break
		if not header_done:
			if c == "\n":
				header_done=True
		else:
			if c != "\n":
				if last_nuc+c == "TA" or last_nuc+c == "ta":
					print(pos, file=out_handle)
				pos+=1
				last_nuc=c

exit()		
				
