import pymysql
import csv

connection = pymysql.connect(host='localhost',
                             user='root',
                             password='root',
                             database='aligner')

with open('examples/blast_results.csv', 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        print(row)
        with connection.cursor() as cursor:
            cursor.execute(
                f"""INSERT INTO blast_cmp (query_sequence, target_sequence, blast_p_value)
                        VALUES ('{row[0]}', '{row[1]}', {row[2].replace(',', '.')})""")
            connection.commit()
