#@title Mount Google Drive
save_to_google_drive = True

from google.colab import drive
drive.mount('/content/drive')

#@markdown ### GoogleAuth
from google.colab import files # Needed in the last cell for packaging and downloading files

from pydrive.drive import GoogleDrive
from pydrive.auth import GoogleAuth
from google.colab import auth
from oauth2client.client import GoogleCredentials
auth.authenticate_user()
gauth = GoogleAuth()
gauth.credentials = GoogleCredentials.get_application_default()
drive = GoogleDrive(gauth)

print("You are logged into Google Drive and are good to go!")


def get_drive():
  # Access Google Drive API for read, write, delete
  service = auth('https://www.googleapis.com/auth/drive.metadata.readonly')

  # Call the Drive v3 API
  results = service.files().list(
    pageSize=10, fields="nextPageToken, files(id, name)").execute()
  items = results.get('files', [])

  if not items:
    print('No files found.')
  else:
    print('Files:')
    for item in items:
      print(u'{0} ({1})'.format(item['name'], item['id']))
