import smtplib
from future.utils import string_types
from email.mime.text import MIMEText

from fanc import config


def send_email(to_address, message, subject='', from_address=None, server=None, credentials=None,
               include_default_subject=True, html=True):
    if subject is None:
        subject = ''

    if include_default_subject:
        subject = '[FAN-C] ' + subject

    if from_address is None:
        if config.email_from_address is None:
            raise ValueError("No from_address. "
                             "Must directly provide email sender information "
                             "or modify fanc.conf file accordingly!")
        from_address = config.email_from_address

    port = None
    if server is None:
        if config.email_smtp_server is None:
            raise ValueError("No server address. "
                             "Must directly provide server address "
                             "or modify fanc.conf file accordingly!")
        server = config.email_smtp_server
        port = config.email_smtp_port
    else:
        if isinstance(server, string_types):
            server = server.split(":")
        server = server[0]
        if len(server) == 2:
            port = server[1]

    if credentials is None or credentials[0] is None or credentials[1] is None:
        username = config.email_smtp_username
        password = config.email_smtp_password
    else:
        username, password = credentials

    # start connection
    if port is None:
        smtp_server = smtplib.SMTP(server)
    else:
        smtp_server = smtplib.SMTP(server, port=port)

    smtp_server.ehlo()
    try:
        smtp_server.starttls()
    except smtplib.SMTPException:
        pass

    if username is not None and password is not None:
        smtp_server.login(username, password)

    if html:
        msg = MIMEText(message, 'html')
    else:
        msg = MIMEText(message, 'plain')

    msg['From'] = from_address
    msg['To'] = to_address
    msg['Subject'] = subject

    smtp_server.sendmail(msg['From'], msg['To'], msg.as_string())
