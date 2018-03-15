from wopmetabarcoding.wrapper.functions import insert_table
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


def create_fastadb(session, model, csv_file, reverse_csv_file, fasta_file):
	"""
	Function creating a fasta file which will be used as a database by vsearch
	:param session: Current session of the database
	:param model: Model of the FileInformation table
	:param csv_file:
	:param reverse_csv_file:
	:return:
	"""
	output_file = open(csv_file, 'w')
	reverse_output_file = open(reverse_csv_file, 'w')
	for line in session.query(model).all():
		if line.tag_forward != "" and line.primer_forward != "" and line.tag_reverse != "" and line.primer_reverse != "":
			if line.file_name == fasta_file:
				output_file.write(">" + line.tag_forward + line.primer_forward)
				output_file.write("\n")
				output_file.write(line.tag_forward + line.primer_forward)
				output_file.write("\n")
			reverse_output_file.write(">" + line.tag_reverse + line.primer_reverse)
			reverse_output_file.write("\n")
			reverse_output_file.write(line.tag_reverse + line.primer_reverse)
			reverse_output_file.write("\n")
	output_file.close()


def read_counter(session, file, model):
	"""
	Function counting occurences of a read in the fasta file and store it in a table
	:param session: Current of the database
	:param file: fasta containing the reads
	:param model: Model of the ReadCount table
	:return: void
	"""
	with open(file, "r") as fasta_file:
		next(fasta_file)
		sequence = ""
		liste_tmp = []
		for line in fasta_file:
			if ">" in line:
				if sequence in liste_tmp:
					session.query(model).filter(model.sequence == sequence).update({model.count: model.count+1})
					sequence = ""
				else:
					obj_readcount = {"sequence": sequence, "count": 1}
					insert_table(session, model, obj_readcount)
					liste_tmp.append(sequence)
					sequence = ""
			else:
				sequence += line.strip()
		if session.query(model.sequence.contains(sequence)) is True and sequence != "":
			session.query(model).filter(model.sequence == sequence).update({model.count: model.count + 1})
			sequence = ""
		else:
			obj_readcount = {"sequence": sequence, "count": 1}
			insert_table(session, model, obj_readcount)
			sequence = ""


def dereplicate(outputcsv, tsv_file):
	"""
	Function use to insert all the data of the vsearch alignment in the database table
	:param session: Current session of the database
	:param model: Model of the table
	:param outputcsv: csv containing the results of the alignment
	:return: void
	"""
	tsv_reads = open(tsv_file, 'w')
	with open(outputcsv, 'r') as file_csv:
		next(file_csv)
		for line in file_csv:
			if line.split("\t")[5] == "1":
				sequence_id = line.split('\t')[0]
				target = line.split('\t')[1]
				tl = line.split('\t')[2]
				qilo = line.split('\t')[3]
				qihi =line.split('\t')[4]
				tilo = line.split('\t')[5]
				tihi = line.split('\t')[6]
				qrow = line.split('\t')[7].strip()
				tag_sequence = ""
				for character in target:
					if character.islower():
						tag_sequence += character.upper()
				if tilo == "1" and tihi == tl and tag_sequence in qrow:
					tsv_reads.write(line)


def fasta_writer(session, model, file_name):
	"""
	Function
	:param session:
	:param model:
	:param file_name:
	:return:
	"""
	name = file_name.replace(".fasta", "_count.fasta")
	file = open(name, "w")
	j = 1
	line = []
	for element in session.query(model).all():
		file.write("> Read number " + str(j) + " count: " + str(element.count))
		file.write("\n")
		# for i in range(len(element.sequence)):
		# 	if (i%80) == 0:
		# 		file.write(sequence)
		# 		file.write('\n')
		file.write(element.sequence + '\n')
		j += 1


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
					sequence = ""
					line = line.replace('>', '')
					id = line.strip().split(" ")[0]
			else:
				line = line.strip()
				sequence += line
			i += 1
		sequence_dico[id] = sequence
		sequence = ""
	return sequence_dico


def insert_read(csv_file, fasta_file, session, file_model):
	filename = fasta_file.replace('.fasta', '_forward_trimmed.csv')
	test_file = open(filename, 'w')
	session.query(file_model).filter(file_model.file_name == fasta_file).update({file_model.forward_trimmed_file: filename})
	reads = read_catcher(fasta_file)
	with open(csv_file, 'r') as csv_file:
		for line in csv_file:
			line_info = line.strip().split('\t')
			read = reads.get(line_info[0])
			qihi = line_info[4]
			trimmed_part = read[0:int(qihi)]
			trimmed_read = read.replace(trimmed_part, "")
			my_read = Seq(trimmed_read, IUPAC.ambiguous_dna)
			reverse_read = my_read.reverse_complement()
			if read is not None:
				line = line.strip()
				test_file.write(line + '\t' + str(reverse_read) + '\n')
			else:
				print(line_info[0])
	test_file.close()


def create_fasta(forward_trimmed_fasta, merged_fasta):
	new_fasta = forward_trimmed_fasta.replace('.csv', '.fasta')
	print(new_fasta)
	csv_file = open(forward_trimmed_fasta, 'r')
	with open(new_fasta, 'w') as fasta_file:
		for line in csv_file:
			line = line.strip()
			line = line.split("\t")
			print(line)
			tag = "".join([character for character in line[1] if character.islower()])
			marker = "".join([character for character in line[1] if character.isupper()])
			fasta_file.write(">" + line[0] + "|" + tag + "|" + marker + "|" + merged_fasta + "\n")
			fasta_file.write(line[8])
			fasta_file.write("\n")
	csv_file.close()


def read_count_variants(session, model1, tsv_file):
	with open(tsv_file, 'r') as file:
		i = 1
		marker = ""
		liste_tmp =[]
		next(file)
		for line in file:
			line_info = line.split('\t')
			if line_info[5] in liste_tmp:
				session.query(model1).filter(model1.read == line_info[5]).update({model1.read_count: model1.read_count + 1})
				sample = '|' + line_info[4]
				data = session.query(model1).filter(model1.read == line_info[5]).first()
				print(data.sample_count)

			else:
				if line_info[2] != marker:
					i = 1
				variant_id = line_info[2] + "_" + line_info[1] + "_" + str(i)
				sample_count = line_info[3] + "-" + line_info[4]
				read_count = 1
				read = line_info[5]
				liste_tmp.append(read)
				obj_obifasta = {'variant_id': variant_id, 'sample_count': sample_count, 'read_count': read_count, 'read': read}
				insert_table(session, model1, obj_obifasta)
				marker = line_info[2]
				i += 1


def attribute_combination(session, fileinformation, model, model2, fasta_file):
	filename = fasta_file.replace('.fasta', '_tmp.tsv')
	file_tmp = open(filename, 'w')
	file_tmp.write("Id" + "\t" + "Marker" + "\t" + "Run" + "\t" + "Sample" + "\t" + "Replicate" + "\t" + "Read" + "\n")
	with open(fasta_file, 'r') as file:
		for line in file:
			if ">" in line:
				line_info = line.strip().split('|')
				id = line_info[0].replace('>', '')
				filename_data = session.query(model).filter(model.sequence_id == id).first()
				data_line = session.query(fileinformation).filter(fileinformation.tag_forward == line_info[1]).filter(fileinformation.primer_forward == line_info[2]).filter(fileinformation.tag_reverse == line_info[3]).filter(fileinformation.primer_reverse == line_info[4]).filter(fileinformation.file_name == filename_data.filename).first()
				file_tmp.write(
					id + "\t" + data_line.run_name + "\t" + data_line.marker_name + "\t" +
					data_line.sample_name + "\t" + data_line.replicate_name + '\t'
				)
			else:
				file_tmp.write(line)
		read_count_variants(session, model2, filename)


def singleton_removing(session, model):
	session.query(model).filter(model.read_count == 1).delete()


def variant_table(session, obifasta, variant):
	for element in session.query(obifasta).all():
		variant_id = element.variant_id
		marker = variant_id.split("_")[0]
		obj_variant = {'variant_id': variant_id, 'marker': marker, 'sequence': element.read}
		insert_table(session, variant, obj_variant)







