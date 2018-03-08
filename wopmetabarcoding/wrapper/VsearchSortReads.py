from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmetabarcoding.wrapper.functions import insert_table, fasta_writer, insert_read, trim_reads, reverse_complement, create_forward_fasta
from sqlalchemy import update
import subprocess


class VsearchSortReads(ToolWrapper):
	__mapper_args__ = {
		"polymorphic_identity": "wopmetabarcoding.wrapper.VsearchSortReads"
	}
	__input_table_fileinformation = "FileInformation"
	__input_table_file = "File"
	__output_file_csvsearch = "csvsearch"
	__output_file_fastadb = "fastadb"
	__output_table_readcount = "ReadCount"
	__output_table_outputcsv = "OutputCsv"
	__output_file_reversefastadb = "reversefastadb"
	__output_file_fowardtrimmed = "fowardtrimmed"

	def specify_input_table(self):
		return [
			VsearchSortReads.__input_table_fileinformation,
			VsearchSortReads.__input_table_file
		]

	def specify_output_file(self):
		return [
			VsearchSortReads.__output_file_csvsearch,
			VsearchSortReads.__output_file_fastadb,
			VsearchSortReads.__output_file_reversefastadb,
			VsearchSortReads.__output_file_fowardtrimmed
		]

	def specify_output_table(self):
		return [
			VsearchSortReads.__output_table_readcount,
			VsearchSortReads.__output_table_outputcsv
		]

	def specify_params(self):
		"""

		:return:
		"""
		return {
			"min_id": "float",
			"minseqlength": "int"
		}

	def create_fastadb(self, session, model, csv_file, reverse_csv_file):
		"""

		:param session:
		:param model:
		:param file:
		:return:
		"""
		output_file = open(csv_file, 'w')
		reverse_output_file = open(reverse_csv_file, 'w')
		for line in session.query(model).all():
			output_file.write(">" + line.tag_forward + line.primer_forward)
			output_file.write("\n")
			output_file.write(line.tag_forward + line.primer_forward)
			output_file.write("\n")
			reverse_output_file.write(">" + line.tag_reverse + line.primer_reverse)
			reverse_output_file.write("\n")
			reverse_output_file.write(">" + line.tag_reverse + line.primer_reverse)
		output_file.close()

	def vsearch_sr(self, file, fasta_db, output_csv):
		"""

		:param file:
		:param fasta_db:
		:param output_csv:
		:return:
		"""
		subprocess.call(
			"vsearch --usearch_global " + file + " --db " + fasta_db + " --id " + str(self.option("min_id")) +
			" --maxhits 1 --maxrejects 0 --maxaccepts 0 --minseqlength " + str(self.option("minseqlength")) +
			" --userout " + output_csv + " --userfields query+target+tl+qilo+qihi+tilo+tihi+qrow", shell=True
		)

	@staticmethod
	def read_counter(session, file, model):
		"""

		:param session:
		:param file:
		:param model:
		:return:
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

	@staticmethod
	def dereplicate(session, model, outputcsv):
		"""

		:param session:
		:param model:
		:param outputcsv:
		:return:
		"""
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
					obj_outputcsv = {
						'sequence_id': sequence_id, 'tl': tl, 'target': target, 'qilo': qilo, 'qihi': qihi, 'tilo': tilo, 'tihi': tihi, 'qrow': qrow
					}
					insert_table(session, model, obj_outputcsv)

	def algo_trim(self, session, output_csv, csv_file, database, fasta_file):
		self.vsearch_sr(fasta_file, database, csv_file)
		self.dereplicate(session, output_csv, csv_file)
		insert_read(session, output_csv, fasta_file)
		trim_reads(session, output_csv)
		reverse_complement(session, output_csv)

	def run(self):
		session = self.session()
		engine = session._WopMarsSession__session.bind
		conn = engine.connect()
		# Output file
		csv_output = self.output_file(VsearchSortReads.__output_file_csvsearch)
		fasta_db = self.output_file(VsearchSortReads.__output_file_fastadb)
		reverse_fasta_db = self.output_file(VsearchSortReads.__output_file_reversefastadb)
		foward_trimmed = self.output_file(VsearchSortReads.__output_file_fowardtrimmed)
		# Table model
		file_information_model = self.input_table(VsearchSortReads.__input_table_fileinformation)
		file_model = self.input_table(VsearchSortReads.__input_table_file)
		readcount_model = self.output_table(VsearchSortReads.__output_table_readcount)
		outputcsv_model = self.output_table(VsearchSortReads.__output_table_outputcsv)
		for element in session.query(file_model).all():
			if element.dereplicate_status is "1":
				self.read_counter(session, element.file_name, readcount_model)
				fasta_writer(session, readcount_model, element.file_name)
			else:
				self.create_fastadb(session, file_information_model, fasta_db, reverse_fasta_db)
				self.algo_trim(session, outputcsv_model, csv_output, fasta_db, element.file_name)
				print("ok")
				create_forward_fasta(session, outputcsv_model, foward_trimmed)
				self.algo_trim(session, outputcsv_model, csv_output, reverse_fasta_db, foward_trimmed)

		session.commit()


