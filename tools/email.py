import os
import smtplib

from email.message import EmailMessage

new_line = '\n'

def send_email(subject, body):
    """Send a simple email."""

    msg = EmailMessage()
    msg.set_content(body)

    msg['Subject'] = subject
    msg['From'] = "cgg-cancer@gu.se" # TODO Get from config
    msg['To'] = "cgg-cancer@gu.se" # TODO Get from config and have different recipients for errors and success
    msg['Cc'] = "hanna.soderstrom@gu.se" # TODO Get from config


    # Send the message
    s = smtplib.SMTP('smtp.gu.se')
    s.send_message(msg)
    s.quit()


def start_email(run_name, samples):
    """Send an email about starting wgs-somatic for samples in a run"""

    subject = f'WGS Somatic start mail {run_name}'

    body = '\n'.join([f'Starting wgs_somatic for the following samples in run {run_name}:',
                      '',
                      f'{new_line}{new_line.join(samples)}',
                      '',
                      'You will get an email when the results are ready.',
                      '',
                      'Best regards,',
                      'CGG Cancer'])

    send_email(subject, body)

def end_email(run_name, samples):
    """Send an email that wgs-somatic has finished running for samples in a run"""

    subject = f'WGS Somatic end mail {run_name}'

    body = '\n'.join([f'WGS somatic has finished successfully for the following samples in run {run_name}:',
                      '',
                      f'{new_line}{new_line.join(samples)}',
                      '',
                      'Best regards,',
                      'CGG Cancer'])

    send_email(subject, body)


#def error_email():




