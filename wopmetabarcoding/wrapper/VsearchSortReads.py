from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmetabarcoding.wrapper.functions import insert_table, fasta_writer
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

	def specify_input_table(self):
		return [
			VsearchSortReads.__input_table_fileinformation,
			VsearchSortReads.__input_table_file
		]

	def specify_output_file(self):
		return [
			VsearchSortReads.__output_file_csvsearch,
			VsearchSortReads.__output_file_fastadb
		]

	def specify_output_table(self):
		return [
			VsearchSortReads.__output_table_readcount
		]

	def specify_params(self):
		return {
			"min_id": "float",
			"minseqlength": "int"
		}

	def create_fastadb(self, session, model, file):
		output_file = open(file, 'w')
		i = 0
		with open(file, 'r') as file_content:
			for line in session.query(model).all():
				if i == 0:
					i += 1
					continue
				output_file.write(">" + line.tag_forward + line.primer_forward)
				output_file.write("\n")
				output_file.write(line.tag_forward + line.primer_forward)
				output_file.write("\n")
		output_file.close()

	def vsearch_sr(self, file, fasta_db, output_csv):
		subprocess.call(
			"vsearch --usearch_global " + file + " --db " + fasta_db + " --id " + str(self.option("min_id")) +
			" --maxhits 1 --maxrejects 0 --maxaccepts 0 --minseqlength " + str(self.option("minseqlength")) +
			" --userout " + output_csv + " --userfields query+target+tl+qilo+qihi+tilo+tihi+qrow", shell=True
		)

	def read_counter(self, session, file, model):
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

	def run(self):
		session = self.session()
		engine = session._WopMarsSession__session.bind
		conn = engine.connect()
		# Output file
		csv_output = self.output_file(VsearchSortReads.__output_file_csvsearch)
		fasta_db = self.output_file(VsearchSortReads.__output_file_fastadb)
		# Table model
		file_information_model = self.input_table(VsearchSortReads.__input_table_fileinformation)
		file_model = self.input_table(VsearchSortReads.__input_table_file)
		readcount_model = self.output_table(VsearchSortReads.__output_table_readcount)
		for element in session.query(file_model).all():
			if element.dereplicate_status is "1":
				self.read_counter(session, element.file_name, readcount_model)
				fasta_writer(session, readcount_model, element.file_name)
			else:
				self.create_fastadb(session, file_information_model, fasta_db)
				# self.vsearch_sr(element.file_name, fasta_db, csv_output)
		session.commit()


