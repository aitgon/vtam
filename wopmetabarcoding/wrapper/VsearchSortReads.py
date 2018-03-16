from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmetabarcoding.wrapper.VsearchSortReads_functions import \
	create_fastadb, dereplicate, read_counter, fasta_writer, insert_read, \
	create_fasta, attribute_combination

import subprocess



class VsearchSortReads(ToolWrapper):
	__mapper_args__ = {
		"polymorphic_identity": "wopmetabarcoding.wrapper.VsearchSortReads"
	}
	# Input
	# Input table
	__input_table_fileinformation = "FileInformation"
	# Input file
	__input_table_file = "File"
	# Output
	# Output file
	__output_file_csvsearch = "csvsearch"
	__output_file_fastadb = "fastadb"
	__output_file_tsv_tmp = "tsv"
	# Output table
	__output_table_readcount = "ReadCount"
	__output_table_obifasta = 'ObiFasta'
	__output_table_variant = 'Variant'

	def specify_input_table(self):
		return [
			VsearchSortReads.__input_table_fileinformation,
			VsearchSortReads.__input_table_file
		]

	def specify_output_file(self):
		return [
			VsearchSortReads.__output_file_csvsearch,
			VsearchSortReads.__output_file_fastadb,
			VsearchSortReads.__output_file_tsv_tmp
		]

	def specify_output_table(self):
		return [
			VsearchSortReads.__output_table_readcount,
			VsearchSortReads.__output_table_obifasta,
			VsearchSortReads.__output_table_variant
		]

	def specify_params(self):
		"""

		:return:
		"""
		return {
			"min_id": "float",
			"minseqlength": "int"
		}

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

	def algo_trim(self, csv_file, database, fasta_file, tsv_file, session, file_model, strain):
		"""

		:param session:
		:param output_csv:
		:param csv_file:
		:param database:
		:param fasta_file:
		:return:
		"""
		self.vsearch_sr(fasta_file, database, csv_file)
		dereplicate(csv_file, tsv_file)
		insert_read(tsv_file, fasta_file, session, file_model, strain)
		session.commit()

	def run(self):
		session = self.session()
		engine = session._WopMarsSession__session.bind
		conn = engine.connect()

		# Input tables models
		file_information_model = self.input_table(VsearchSortReads.__input_table_fileinformation)
		file_model = self.input_table(VsearchSortReads.__input_table_file)

		# Output file
		csv_output = self.output_file(VsearchSortReads.__output_file_csvsearch)
		fasta_db = self.output_file(VsearchSortReads.__output_file_fastadb)
		tsv_tmp = self.output_file(VsearchSortReads.__output_file_tsv_tmp)

		# Output tables models
		readcount_model = self.output_table(VsearchSortReads.__output_table_readcount)
		obifasta_model = self.output_table(VsearchSortReads.__output_table_obifasta)
		variant_model = self.output_table(VsearchSortReads.__output_table_variant)

		for element in session.query(file_model).all():
			if element.dereplicate_status is "1":
				read_counter(session, element.file_name, readcount_model)
				fasta_writer(session, readcount_model, element.file_name)
			else:
				create_fastadb(session, file_information_model, fasta_db, element.file_name)
				self.algo_trim(csv_output, fasta_db, element.file_name, tsv_tmp, session, file_model, 'forward')
				create_fasta(element.forward_trimmed_file)
		session.commit()
		for element in session.query(file_model).all():
			forward_trimmed_fasta = element.forward_trimmed_file.replace('.csv', '.fasta')
			create_fastadb(session, file_information_model, fasta_db, forward_trimmed_fasta)
			self.algo_trim(csv_output, fasta_db, forward_trimmed_fasta, tsv_tmp, session, file_model, 'reverse')
			create_fasta(element.output_reverse_file)
		session.commit()
		for element in session.query(file_model).all():
			attribute_combination(session, file_information_model, file_model, element.output_reverse_file, element.file_name)
		# for element in session.query(file_model).all():
