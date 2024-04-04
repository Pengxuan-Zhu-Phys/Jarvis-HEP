#!/usr/bin/env python3 

import sqlite3
import json
import csv 
from base import Base
from pprint import pprint
import pandas as pd

class SampleDatabase(Base):
    def __init__(self, db_path):
        self.db_path = db_path
        self._initialize_database()

    def _initialize_database(self):
        """Initialize the database with the necessary tables."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        # Create table if it doesn't exist
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS samples (
                uuid TEXT PRIMARY KEY,
                data TEXT NOT NULL
            )
        ''')
        conn.commit()
        conn.close()

    def insert_sample_row(self, sample_info):
        """Insert a sample into the database."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        try:
            # print("Line 32 ->")
            # pprint(sample_info)
            # Serialize observables to JSON for storage
            data_json = json.dumps(sample_info)
            cursor.execute('''
                INSERT INTO samples (uuid, data) VALUES (?, ?)
            ''', (sample_info["uuid"], data_json))
            conn.commit()
        except sqlite3.IntegrityError:
            print(f"Sample with UUID {sample_info['uuid']} already exists in the database.")
        finally:
            conn.close()

    def get_sample_by_uuid(self, uuid):
        """Retrieve a sample from the database by its UUID."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute('''
            SELECT data FROM samples WHERE uuid=?
        ''', (uuid,))
        row = cursor.fetchone()
        conn.close()
        if row:
            # Deserialize data from JSON
            data = json.loads(row[0])
            return data
        else:
            return None

    def export_database_to_csv(self, csv_path):
        """Export the database data to a CSV file."""
        conn = sqlite3.connect(self.db_path)
        # cursor = conn.cursor()

        # # Fetch all records from the database
        # cursor.execute('SELECT uuid, data FROM samples')
        # records = cursor.fetchall()

        # # Open the CSV file for writing
        # with open(csv_path, 'w', newline='') as csv_file:
        #     writer = csv.writer(csv_file)
        #     # Write the header row
        #     writer.writerow(['UUID', 'Data'])

        #     for uuid, data_json in records:
        #         # Deserialize the JSON data
        #         data = json.loads(data_json)
        #         # Flatten the data into a single row, if necessary
        #         # For simplicity, we'll just write the UUID and the serialized JSON
        #         writer.writerow([uuid, json.dumps(data)])

        # conn.close()


    # Execute the query to fetch all records from the database
        df = pd.read_sql_query('SELECT uuid, data FROM samples', conn)

        # Close the connection to the database
        conn.close()

        # Normalize and flatten the JSON data in the 'data' column
        # Note: This assumes that the 'data' column contains serialized JSON strings
        records = df['data'].apply(json.loads).tolist()  # Deserialize JSON strings to dicts
        flattened_data = pd.json_normalize(records)  # Flatten the JSON structure

        # Include the 'uuid' column from the original DataFrame
        flattened_data['UUID'] = df['uuid']

        # Reorder the columns to have 'UUID' as the first column
        cols = ['UUID'] + [col for col in flattened_data.columns if col != 'uuid']
        flattened_data = flattened_data[cols]

        # Export the flattened DataFrame to CSV
        flattened_data.to_csv(csv_path, index=False)
            # print(f"Database export completed to {csv_path}")
