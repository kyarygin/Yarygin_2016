from zipfile import ZipFile
import requests
import sys
import os

try:
    input_ = raw_input
except NameError:
    input_ = input


filename2id = {
    'BGI_blastdb.zip': '0Byeru5YHfXc-R1lPWVlPd19vYnM',
    'BGI_coverage.zip': '0Byeru5YHfXc-QVJFNWZycEJmWEE',
    'blast_results.zip': '0Byeru5YHfXc-ME5KbVYwVlNzX1k',
    'gene_groups_abund.zip': '0Byeru5YHfXc-enlvU1pycmJ4c1E',
    'index.zip': '0Byeru5YHfXc-OTdteGUxU05HSlU',
    'metadata.zip': '0Byeru5YHfXc-NFJQRVItUWFJU1k'
}

filename2size = {
    'BGI_blastdb.zip': 694118894,
    'BGI_coverage.zip': 1040855156,
    'blast_results.zip': 427363616,
    'gene_groups_abund.zip': 50421398,
    'index.zip': 8628492365,
    'metadata.zip': 34715942
}

def human_readable(size):
    for unit in ['', 'K', 'M', 'G', 'T']:
        if abs(size) < 1024.:
            return '{:3.1f}{}'.format(size, unit)
        size /= 1024.
    return '{:.1f}{}'.format(size, unit)

def download_file_from_google_drive(id, filename):
    url = 'https://docs.google.com/uc?export=download'

    session = requests.Session()

    response = session.get(url, params={'id': id}, stream=True)
    token = get_confirm_token(response)

    if token:
        params = {'id': id, 'confirm': token}
        response = session.get(url, params=params, stream=True)

    save_response_content(response, filename)

def get_confirm_token(response):
    for key, value in response.cookies.items():
        if key.startswith('download_warning'):
            return value
    return None

def save_response_content(response, filename):
    CHUNK_SIZE = 1024 * 32
    downloaded_size = 0
    with open(filename, 'wb') as f:
        for chunk in response.iter_content(CHUNK_SIZE):
            if chunk:
                f.write(chunk)
                downloaded_size += CHUNK_SIZE
                msg = 'Downloading {}: {} / {}     \r'.format(filename,
                    human_readable(downloaded_size),
                    human_readable(filename2size[filename])
                )
                sys.stdout.write(msg)
                sys.stdout.flush()
    sys.stdout.write('\n')

def zip_extract(filename):
    sys.stdout.write('Extracting {} ... '.format(filename))
    sys.stdout.flush()
    z_file = ZipFile(filename)
    z_file.extractall()
    sys.stdout.write('Done\n')
    os.remove(filename)

if __name__ == '__main__':
    sys.stdout.write('This script will download {} on your computer, continue? [Y/n]\n'.format(human_readable(sum(filename2size.values()))))
    resp = input_()
    if resp not in set(['Y', 'y', 'Yes', 'yes', '']):
        sys.exit()
    for filename, file_id in filename2id.items():
        download_file_from_google_drive(file_id, filename)
        zip_extract(filename)

