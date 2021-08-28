# @title Procedure time_check
"""
Prints current time and remaining disk as a progress check
"""
from datetime import datetime
import shutil

from dtumemutils import track

@track
def time_check(message_str):
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print(message_str + current_time)

    total, used, free = shutil.disk_usage("/")
    print("Free: %d GiB" % (free // (2 ** 30)))
