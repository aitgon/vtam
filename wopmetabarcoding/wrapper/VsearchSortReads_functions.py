from wopmetabarcoding.wrapper.functions import insert_table
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import subprocess
import sqlite3


def create_fastadb(session, model, csv_file, fasta_file):
	"""
	Function creating a fasta file which will be used as a database by vsearch
	:param session: Current session of the database
	:param model: Model of the FileInformation table
	:param csv_file:
	:param reverse_csv_file:
	:return:
	"""
	output_file = open(csv_file, 'w')
	for line in session.query(model).all():
		if line.tag_forward != "" and line.primer_forward != "" and line.tag_reverse != "" and line.primer_reverse != "":
			if line.file_name == fasta_file:
				output_file.write(">" + line.tag_forward + line.primer_forward)
				output_file.write("\n")
				output_file.write(line.tag_forward + line.primer_forward)
				output_file.write("\n")
			else:
				output_file.write(">" + line.tag_reverse + line.primer_reverse)
				output_file.write("\n")
				output_file.write(line.tag_reverse + line.primer_reverse)
				output_file.write("\n")
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
		file.write(element.sequence + '\n')
		j += 1


def read_catcher(conn, fasta_file):
	try:
		conn.execute("DROP TABLE IF EXISTS reads_fasta")
		conn.execute("CREATE TABLE  reads_fasta (id VARCHAR, seq VARCHAR)")
		for record in SeqIO.parse(fasta_file, 'fasta'):
			conn.execute("INSERT INTO reads_fasta (id, seq) VALUES (?, ?)", (str(record.description.split()[0]), str(record.seq)))
		conn.commit()
	except UnicodeDecodeError:
		pass


def insert_read(csv_file, fasta_file, session, file_model, strain):
	if strain == 'forward':
		print(fasta_file)
		filename = fasta_file.replace('.fasta', '_forward_trimmed.csv')
		session.query(file_model).filter(file_model.file_name == fasta_file).update({file_model.forward_trimmed_file: filename})
	else:
		fasta_csv = fasta_file.replace('.fasta', '.csv')
		filename = fasta_file.replace('.fasta', '_output_reverse_trimmed.csv')
		session.query(file_model).filter(file_model.forward_trimmed_file == fasta_csv).update({file_model.output_reverse_file: filename})
	test_file = open(filename, 'w')
	conn = sqlite3.connect('db.sqlite')
	read_catcher(conn, fasta_file)
	with open(csv_file, 'r') as csv_file:
		for line in csv_file:
			line_info = line.strip().split('\t')
			read_cursor = conn.execute('SELECT seq FROM reads_fasta WHERE id=?', (line_info[0],))
			for row in read_cursor:
				read = row[0]
			read_cursor.close()
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
	conn.close()


def create_fasta(forward_trimmed_fasta):
	new_fasta = forward_trimmed_fasta.replace('.csv', '.fasta')
	print(new_fasta)
	csv_file = open(forward_trimmed_fasta, 'r')
	with open(new_fasta, 'w') as fasta_file:
		for line in csv_file:
			line = line.strip()
			line = line.split("\t")
			tag = "".join([character for character in line[1] if character.islower()])
			marker = "".join([character for character in line[1] if character.isupper()])
			fasta_file.write(">" + line[0] + "|" + tag + "|" + marker + "\n")
			fasta_file.write(line[8])
			fasta_file.write("\n")
	csv_file.close()


def attribute_combination(session, model, model2, csv_file, filename):
	output_filename = csv_file.replace('_forward_trimmed_output_reverse_trimmed.csv', '_combination.tsv')
	session.query(model2).filter(model2.file_name == filename).update({model.final_csv: output_filename})
	output = open(output_filename, 'w')
	with open(csv_file, 'r') as csv_input:
		for line in csv_input:
			line = line.strip()
			line = line.split('\t')
			id = line[0]
			id = id.strip()
			id = id.split('|')
			tag_reverse = "".join([character for character in line[1] if character.islower()])
			data = session.query(model).filter(model.tag_forward == id[1]).filter(model.tag_reverse == tag_reverse).filter(model.file_name == filename).first()
			output.write(
				id[0] + "\t" + data.marker_name+ "\t" +data.tag_forward + "\t" + data.tag_reverse + "\t" + data.sample_name
				+ "\t" + data.replicate_name + "\t" + line[8] + "\n"
			)