import sys
import csv


spamWriter = csv.writer(open('eggs.csv', 'wb'), delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
spamWriter.writerow(['Spam'] * 5 + ['Baked Beans'])
spamWriter.writerow(['Spam', 'Lovely Spam', 'Wonderful Spam'])
