from wopmetabarcoding.wrapper.functions import insert_table


def insert_marker(session, model, line):
    """
    Function parsing the line to obtain element to insert in Marker table
    :param session: Current session of the database
    :param model: Model of the Marker table
    :param line: Line of the csv file
    :return: void
    """
    marker_name = line.split(',')[4]
    obj_marker = {'marker_name': marker_name}
    insert_table(session, model, obj_marker)


def insert_primer(session, model, line):
    """
    Function parsing the line to obtain element to insert in PrimerPair table
    :param session: Current session of the database
    :param model: Model of the PrimerPair table
    :param line: Line of the csv file
    :return: void
    """
    primer_forward = line.split(',')[1]
    primer_reverse = line.split(',')[3]
    obj_primer = {'primer_forward': primer_forward, 'primer_reverse': primer_reverse}
    insert_table(session, model, obj_primer)


def insert_tagpair(session, model, line):
    """
    Function parsing the line to obtain element to insert in TagPair table
    :param session: Current session of the database
    :param model: Model of the TagPair table
    :param line: Line of the csv file
    :return: void
    """
    tag_forward = line.split(',')[0]
    tag_reverse = line.split(',')[2]
    obj_tag = {'tag_forward': tag_forward, 'tag_reverse': tag_reverse}
    insert_table(session, model, obj_tag)


def insert_file(session, model, line):
    """
    Function parsing the line to obtain element to insert in File table
    :param session: Current session of the database
    :param model: Model of the File table
    :param line: Line of the csv file
    :return: void
    """
    file_name = line.split(',')[7]
    run_name = line.split(',')[8].strip()
    if line.split(',')[0] == "" or line.split(',')[2] == "" or line.split(',')[3] == "" or line.split(',')[4] == "":
        dereplicate = True
    else:
        dereplicate = False
    obj_file = {'name': file_name, 'run_name': run_name, 'trimmed_status': dereplicate}
    insert_table(session, model, obj_file)


def insert_sample(session, model, line):
    """
    Function parsing the line to obtain element to insert in Biosample table
    :param session: Current session of the database
    :param model: Model of the Biosample table
    :param line: Line of the csv file
    :return: void
    """
    sample_name = line.split(',')[5]
    obj_sample = {'name': sample_name}
    insert_table(session, model, obj_sample)

def insert_replicate(session, model, line):
    biosample_name = line.split(',')[5]
    marker_name = line.split(',')[4]
    file_name = line.split(',')[7]
    replicate_name = line.split(',')[6]
    obj_replicate = {'biosample_name': biosample_name, 'marker_name': marker_name, 'file_name': file_name, 'name': replicate_name}
    insert_table(session, model, obj_replicate)


def insert_replicate_marker(session, model, line):
    replicate_name = line.split(',')[6]
    marker_name = line.split(',')[4]
    obj_replicatemarker = {'name': replicate_name, 'marker_name': marker_name}
    insert_table(session, model, obj_replicatemarker)

def insert_fileinformation(session, model, line):
    """
    Function parsing the line to obtain element to insert in SampleInformation table
    :param session: Current session of the database
    :param model: Model of the SampleInformation table
    :param line: Line of the csv file
    :return: void
    """
    marker_name = line.split(',')[4]
    primer_forward = line.split(',')[1]
    primer_reverse = line.split(',')[3]
    tag_forward = line.split(',')[0]
    tag_reverse = line.split(',')[2]
    file_name = line.split(',')[7]
    run_name = line.split(',')[8].strip()
    sample_name = line.split(',')[5]
    replicate_name = line.split(',')[6]
    obj_fileinformation = {
        'marker_name': marker_name, 'tag_forward': tag_forward, 'tag_reverse': tag_reverse,
        'primer_forward': primer_forward, 'primer_reverse': primer_reverse, 'file_name': file_name,
        'run_name': run_name, 'sample_name': sample_name, 'replicate_name': replicate_name,
    }
    insert_table(session, model, obj_fileinformation)