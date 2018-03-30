def insert_table(session, model, obj):
    try:  # checks if exists Phenotype in db
        session.query(model).filter_by(**obj).one()
    except:  # if not, add
        session.add(model(**obj))

