# @markdown ### GoogleAuth
from pydrive.drive import GoogleDrive
from pydrive.auth import GoogleAuth
from google.colab import auth
from oauth2client.client import GoogleCredentials

def authorise():
    # @markdown ### GoogleAuth
    from pydrive.drive import GoogleDrive
    from pydrive.auth import GoogleAuth
    from google.colab import auth
    from oauth2client.client import GoogleCredentials
    auth.authenticate_user()
    gauth = GoogleAuth()
    gauth.credentials = GoogleCredentials.get_application_default()
    return gauth

def mount_drive(gauth, location: str):
    drive = GoogleDrive(gauth)
    print("You are logged into Google Drive and are good to go!")


def main2():
    gauth= authorise()
    drive = GoogleDrive(gauth)
    print("You are logged into Google Drive and are good to go!")


if __name__ == "__main__":
main2()
