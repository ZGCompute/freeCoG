# PyEmail_Me.py
# Author: Zack Greenberg
# UCSF Neuroimaging Center
# Department of Neurology
# Last Editd May 1st, 2014

# This script can be used to automate emails to yourself and others for important reminders

#!/usr/bin/python
 
# Import smtplib for the actual sending function
import smtplib
 
# For guessing MIME type
import mimetypes
 
# Import the email modules we'll need
import email
import email.mime.application
 
#Import sys to deal with command line arguments
import sys
 
# Create a text/plain message
msg = email.mime.Multipart.MIMEMultipart()
msg['Subject'] = 'Imaging CORE Update'
msg['From'] = '*MY EMAIL*'
recipients = ['*LIST OF EMAILS TO SEND TO*']
msg['To'] = ", ".join(recipients)
 
# The main body is just another attachment
body = email.mime.Text.MIMEText("""All,
 
-Zack
--
Zachary Greenberg
Imaging CORE Manager
UCSF Memory And Aging Center
Department of Neurology
Email: """)
msg.attach(body)
 
# PDF attachment block code
directory='/Users/zackgreenberg/Desktop/MAC/Figures/MAC_Scan_tracking.png' #sys.argv[1] to get current dir
 
# Split de directory into fields separated by / to substract filename
spl_dir=directory.split('/')
 
# We attach the name of the file to filename by taking the last
# position of the fragmented string, which is, indeed, the name
# of the file we've selected 
filename=spl_dir[len(spl_dir)-1]
 
# We'll do the same but this time to extract the file format (pdf, epub, docx...)
spl_type=directory.split('.')
 
type=spl_type[len(spl_type)-1]
 
fp=open(directory,'rb')
att = email.mime.application.MIMEApplication(fp.read(),_subtype=type)
fp.close()
att.add_header('Content-Disposition','attachment',filename=filename)
msg.attach(att)
 
# send via Gmail server
s = smtplib.SMTP('smtp.gmail.com:587')
s.starttls()
s.login('*EMAIL*','*PASSWORD*')
s.sendmail('*EMAIL*',recipients, msg.as_string())
s.quit()
