from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from sqlalchemy.sql import select
import subprocess


class VsearchSortReads(ToolWrapper):
	__mapper_args__ = {
		"polymorphic_identity": "wopmetabarcoding.wrapper.VsearchSortReads"
	}
	__input_table_fileinformation = "FileInformation"
	__input_table_file = "File"
	__output_file_csvsearch = "csvsearch"
	__output_file_fastadb = "fastadb"

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

	def vsearch_sr(self, session, file_model, fasta_db, output_csv):
		for line in session.query(file_model).all():
			subprocess.call(
				"vsearch --usearch_global " + line.file_name + " --db " + fasta_db + " --id " + str(self.option("min_id")) +
				" --maxhits 1 --maxrejects 0 --maxaccepts 0 --minseqlength " + str(self.option("minseqlength")) +
				" --userout " + output_csv + " --userfields query+target+tl+qilo+qihi+tilo+tihi+qrow", shell=True
			)

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
		self.create_fastadb(session, file_information_model, fasta_db)
		self.vsearch_sr(session, file_model, fasta_db, csv_output)


