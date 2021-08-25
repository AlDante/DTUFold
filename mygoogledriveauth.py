from __future__ import print_function
import pickle
import os.path
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request

'''
Google Drive API Scopes: (See https://developers.google.com/identity/protocols/googlescopes#drivev3)
https://www.googleapis.com/auth/drive	                    See, edit, create, and delete all of your Google Drive files
https://www.googleapis.com/auth/drive.appdata	            View and manage its own configuration data in your Google Drive
https://www.googleapis.com/auth/drive.file	                View and manage Google Drive files and folders that you have opened or created with this app
https://www.googleapis.com/auth/drive.metadata	            View and manage metadata of files in your Google Drive
https://www.googleapis.com/auth/drive.metadata.readonly	View metadata for files in your Google Drive
https://www.googleapis.com/auth/drive.photos.readonly	    View the photos, videos and albums in your Google Photos
https://www.googleapis.com/auth/drive.readonly	            See and download all your Google Drive files
https://www.googleapis.com/auth/drive.scripts	            Modify your Google Apps Script scripts' behavior
'''

# If modifying the scope, first delete the file token.pickle.
def auth(scope):

    creds = None
    # The file token.pickle stores the user's access and refresh tokens, and is
    # created automatically when the authorization flow completes for the first
    # time.
    if os.path.exists('token.pickle'):
        with open('token.pickle', 'rb') as token:
            creds = pickle.load(token)
    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(
                'credentials.json', scope)
            creds = flow.run_local_server(port=0)
        # Save the credentials for the next run
        with open('token.pickle', 'wb') as token:
            pickle.dump(creds, token)

    service = build('drive', 'v3', credentials=creds)
    return(service)


if __name__ == '__main__':
    # To authorize full drive access use 'https://www.googleapis.com/auth/drive'
    service = auth('https://www.googleapis.com/auth/drive.metadata.readonly')

    # Call the Drive v3 API
    results = service.files().list(pageSize=10, fields="nextPageToken, files(id, name)").execute()
    items = results.get('files', [])

    if not items:
        print('No files found.')
    else:
        print('Files:')
        for item in items:
            print(u'{0} ({1})'.format(item['name'], item['id']))
