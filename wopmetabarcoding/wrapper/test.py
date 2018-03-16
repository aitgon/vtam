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