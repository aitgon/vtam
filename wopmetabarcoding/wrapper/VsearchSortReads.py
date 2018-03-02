from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from sqlalchemy.sql import select


class VsearchSortReads(ToolWrapper):
	__mapper_args__ = {
		"polymorphic_identity": "wopmetabarcoding.wrapper.VsearchSortReads"
	}
	__input_table_fileinformation = "FileInformation"
	__output_file_csvsearch = "csvsearch"
	__output_file_fastadb = "fastadb"

	def specify_input_table(self):
		return [
			VsearchSortReads.__input_table_fileinformation
		]

	def specify_output_file(self):
		return [
			VsearchSortReads.__output_file_csvsearch,
			VsearchSortReads.__output_file_fastadb
		]

	def create_fastadb(self, session, model, file):
		output_file = open('fasta_db1.fasta', 'w')
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

	def run(self):
		session = self.session()
		engine = session._WopMarsSession__session.bind
		conn = engine.connect()
		# Output file
		csv_output = self.output_file(VsearchSortReads.__output_file_csvsearch)
		fasta_db = self.output_file(VsearchSortReads.__output_file_fastadb)
		# Table model
		file_information_model = self.input_table(VsearchSortReads.__input_table_fileinformation)
		self.create_fastadb(session, file_information_model, fasta_db)


