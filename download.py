import requests

def download_file_from_google_drive(id, destination):
    url = 'https://docs.google.com/uc?export=download'

    session = requests.Session()

    response = session.get(url, params={'id': id}, stream=True)
    token = get_confirm_token(response)

    if token:
        params = {'id': id, 'confirm': token}
        response = session.get(url, params=params, stream=True)

    save_response_content(response, destination)

def get_confirm_token(response):
    for key, value in response.cookies.items():
        if key.startswith('download_warning'):
            return value
    return None

def save_response_content(response, destination):
    CHUNK_SIZE = 1024 * 32
    with open(destination, 'wb') as f:
        for chunk in response.iter_content(CHUNK_SIZE):
            if chunk:
                f.write(chunk)

if __name__ == '__main__':
    file_id = '0Byeru5YHfXc-enlvU1pycmJ4c1E'
    destination = './gene_groups_abund.zip'
    download_file_from_google_drive(file_id, destination)

filename2id = {
    'BGI_blastdb.zip': '0Byeru5YHfXc-R1lPWVlPd19vYnM',
    'BGI_coverage.zip': '0Byeru5YHfXc-QVJFNWZycEJmWEE',
    'blast_results.zip': '0Byeru5YHfXc-ME5KbVYwVlNzX1k',
    'gene_groups_abund.zip': '0Byeru5YHfXc-enlvU1pycmJ4c1E',
    'index.zip': '0Byeru5YHfXc-anQwODlXcU5nTGc'
}




