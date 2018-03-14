from wopmetabarcoding.wrapper.functions import insert_table
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


def create_fastadb(session, model, csv_file, reverse_csv_file):
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
	tsv_reads = open(tsv_file)
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
				if tilo != "1":
					raise Exception("The first position of the alignment is 1")
				if tihi != tl:
					raise Exception("The last position of the alignment didn't represent the length of the target")
				tag_sequence = ""
				for character in target:
					if character.islower():
						tag_sequence += character
				if tag_sequence not in qrow:
					raise Exception("No match between the tag and the alignment")
				tsv_file.write("")

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


def insert_read(session, model, file):
	"""

	:param session:
	:param model:
	:param file:
	:return:
	"""
	with open(file, 'r') as file_fasta:
		i = 1
		sequence = ""
		liste_reads = []
		for line in file_fasta:
			if ">" in line:
				if i == 1:
					sequence_id = line.strip()
					sequence_id = sequence_id.replace('>', '')
					continue
				else:
					liste_reads.append({'sequence_id': sequence_id, 'sequence': sequence})
					sequence = ""
					sequence_id = line.strip()
					sequence_id = sequence_id.replace('>', '')
			else:
				seq = line.replace('\n', '')
				sequence += seq
			i += 1
	liste_reads.append({'sequence_id': sequence_id, 'sequence': sequence})
	for dico in liste_reads:
		for element in session.query(model).all():
			if element.sequence_id in dico.get("sequence_id"):
				session.query(model).filter(model.sequence_id == element.sequence_id).update({model.read: dico.get("sequence"), model.filename: file})


def trim_reads(session, model, file):
	"""
	Function trim the reads of the tag and primer
	:param session:
	:param model:
	:return:
	"""
	for element in session.query(model).all():
		if element.filename == file:
			read_sequence = element.read
			target = element.read[0:element.qihi]
			trimmed_read = read_sequence.replace(target, '')
			session.query(model).filter(model.sequence_id == element.sequence_id).update({model.read: trimmed_read})


def reverse_complement(session, model):
	"""
	Function returning the reverse complement of the read to the database
	:param session: Current session of the database
	:param model: Model of the table
	:return: void
	"""
	for element in session.query(model).all():
		my_read = Seq(element.read, IUPAC.ambiguous_dna)
		reverse_read = my_read.reverse_complement()
		session.query(model).filter(model.sequence_id == element.sequence_id).update({model.reverse: str(reverse_read)})


def create_fasta(session, model, file):
	with open(file, 'w') as fasta_file:
		for line in session.query(model).all():
			tag = "".join([character for character in line.target if character.islower()])
			marker = "".join([character for character in line.target if character.isupper()])
			fasta_file.write(">" + line.sequence_id + "|" + tag + "|" + marker + "\n")
			fasta_file.write(line.reverse)
			fasta_file.write("\n")


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







