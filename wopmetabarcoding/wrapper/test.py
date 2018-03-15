def read_catcher(fasta_file):
	sequence_dico = {}
	sequence = ""
	i = 1
	with open(fasta_file, 'r') as fasta_file:
		for line in fasta_file:
			if ">" in line:
				if i == 1:
					line = line.replace('>', '')
					id = line.strip().split(" ")[0]
				else:
					sequence_dico[id] = sequence
			else:
				line = line.strip()
				sequence += line
			i += 1
	return sequence_dico


def insert_read(csv_file, fasta_file):
	filename = fasta_file.replace('.fasta', '_csv.csv')
	test_file = open(filename, 'w')
	reads = read_catcher(fasta_file)
	with open(csv_file, 'r') as csv_file:
		for line in csv_file:
			line_info = line.strip().split('\t')
			read = reads.get(line_info[0])
			test_file.write(line + '\t' + read + '\n')
