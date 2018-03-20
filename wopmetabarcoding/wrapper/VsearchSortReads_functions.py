from wopmetabarcoding.wrapper.functions import insert_table
import sqlalchemy
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
		conn.execute("CREATE TABLE  reads_fasta (id VARCHAR PRIMARY KEY , seq VARCHAR)")
		for record in SeqIO.parse(fasta_file, 'fasta'):
			conn.execute("INSERT INTO reads_fasta (id, seq) VALUES (?, ?)", (str(record.description.split()[0]), str(record.seq)))
		conn.commit()
	except UnicodeDecodeError:
		pass


def insert_read(csv_file, fasta_file, session, file_model, strain):
	"""

	:param csv_file:
	:param fasta_file:
	:param session:
	:param file_model:
	:param strain:
	:return:
	"""
	if strain == 'forward':
		# Names the files and the database for forward trim purpose
		filename = fasta_file.replace('.fasta', '_forward_trimmed.csv')
		database_name = filename.replace('.csv', '.sqlite')
		session.query(file_model).filter(file_model.file_name == fasta_file).update({file_model.forward_trimmed_file: filename})
	else:
		# Names the files and the database for reverse trim purpose
		fasta_csv = fasta_file.replace('.fasta', '.csv')
		filename = fasta_file.replace('.fasta', '_output_reverse_trimmed.csv')
		database_name = filename.replace('.csv', '.sqlite')
		session.query(file_model).filter(file_model.forward_trimmed_file == fasta_csv).update({file_model.output_reverse_file: filename})
	test_file = open(filename, 'w')
	# Creating a temp database file to enhance speed and stock Id <-> read combination
	conn = sqlite3.connect(database_name)
	# conn.execute("PRAGMA SYNCHRONOUS = OFF ")
	# Function used to insert Id <-> read combination
	read_catcher(conn, fasta_file)
	with open(csv_file, 'r') as csv_file:
		for line in csv_file:
			line_info = line.strip().split('\t')
			read_cursor = conn.execute('SELECT seq FROM reads_fasta WHERE id=?', (line_info[0],))
			read_list = list(read_cursor.fetchone())
			read_cursor.close()
			read = read_list[0]
			qihi = line_info[4]
			trimmed_part = read[0:int(qihi)]
			trimmed_read = read.replace(trimmed_part, "")
			my_read = Seq(trimmed_read, IUPAC.ambiguous_dna)
			reverse_read = my_read.reverse_complement()
			line = line.strip() + '\t' + str(reverse_read) + '\n'
			test_file.write(line)
	test_file.close()


def create_fasta(forward_trimmed_fasta):
	new_fasta = forward_trimmed_fasta.replace('.csv', '.fasta')
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
	session.query(model2).filter(model2.file_name == filename).update({model2.final_csv: output_filename})
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


def count_reads(session, model):
	"""
	Function allowing to count a reads for variants
	:param session: current session of the database
	:param model:
	:return:
	"""
	for element in session.query(model).all():
		conn_name = element.marker_name + ".sqlite"
		print(conn_name)
		session.query(model).filter(model.marker_name == element.marker_name).update({model.db_marker: conn_name})
		conn = sqlite3.connect(conn_name)
		conn.execute("DROP TABLE IF EXISTS count_read")
		conn.execute("CREATE TABLE  count_read (id VARCHAR PRIMARY KEY , marker VARCHAR, count INT, seq VARCHAR)")
		with open(element.marker_file) as marker_file:
			i = 1
			for line in marker_file:
				line = line.strip()
				line = line.split('\t')
				check_read = conn.execute('SELECT EXISTS (SELECT id FROM count_read WHERE seq=?)', (line[6],))
				for row in check_read.fetchone():
					if row != 0:
						conn.execute('UPDATE count_read SET count = count + 1 WHERE seq=?', (line[6],))
					else:
						variant_id = line[1] + "_variant_" + str(i)
						marker_name = line[1]
						read_count = 1
						sequence = line[6]
						conn.execute("INSERT INTO count_read (id, marker, count, seq) VALUES (?, ?, ?, ?)", (variant_id, marker_name, int(read_count), sequence))
						i += 1
		check_read.close()
		conn.execute('DELETE FROM count_read WHERE count=1')
		conn.commit()
		conn.close()
		session.commit()


def gather_files(session, marker_model, file_model):
	for element in session.query(marker_model).all():
		filename = element.marker_name + "_file"
		file_place = "data/" + filename + ".csv"
		session.query(marker_model).filter(marker_model.marker_name == element.marker_name).update({marker_model.marker_file: file_place})
		with open(file_place, 'w') as csv_file:
			for things in session.query(file_model).all():
				if element.marker_name in things.final_csv:
					with open(things.final_csv, 'r') as input_file:
						csv_file.write(input_file.read())
		session.commit()


def insert_variant(session, marker_model, variant_model):
	for element in session.query(marker_model).all():
		conn = sqlite3.connect(element.db_marker)
		variant_data = conn.execute("SELECT id, marker, seq FROM count_read")
		for row in variant_data.fetchall():
			obj_variant = {'variant_id': row[0], 'marker': row[1], 'sequence': row[2]}
			insert_table(session, variant_model, obj_variant)











