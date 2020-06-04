#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python
import os
import json
import time

def read_config():
    with open("run_wrapper_scripts/configs/run_config.json", 'r') as config:
        config_data = json.load(config)
        return config_data

def read_investigators():
    with open("run_wrapper_scripts/configs/investigator_list.json", 'r') as config:
        config_data = json.load(config)
        return config_data


def get_time():
    nowtime = time.strftime("%Y-%m-%d-%H-%M-%S")
    return nowtime

def logger(message, runid, sampleid=False):
    config = read_config()
    log_location = config["start"]["logloc"]
    current_date = time.strftime("%Y-%m-%d")
    log_location = f"{log_location}/{runid}"
    if not os.path.isdir(log_location):
        os.mkdir(log_location)
    if sampleid:
        log_location = f"{log_location}/{runid}/{sampleid}"
        if not os.path.isdir(log_location):
            os.mkdir(log_location)
        logname = f"{log_location}/{sampleid}_{current_date}"
    else:
        logname = f"{log_location}/{runid}_{current_date}"
    if type(message) is dict:
        with open(f"{logname}.json", "w") as jsondump:
            json.dump(message, jsondump, ensure_ascii=False, indent=4)
    else:
        logname = f"{logname}.log"
        logfile = open(logname, "a+")
        logfile.write(f"{get_time()}: {message}" + "\n")
        print(f"{get_time()}: {message}")
        logfile.close()


def send_mail(subj, body, bifogade_filer, recipients=[], sender='clinicalgenomics@gu.se'):
    from marrow.mailer import Mailer, Message
    """Send mail to specified recipients."""
    recipients = [*recipients]
    mailer = Mailer(dict(
        transport=dict(use='smtp',
                       host='smtp.gu.se')))

    mailer.start()
    message = Message(
        subject=f'{subj}',
        plain=f'{body}',
        author=sender,
        to=recipients,)
    if bifogade_filer:
        for fil in bifogade_filer:
            if not os.path.isfile(fil):
                message.plain += "\n"
                message.plain += f"Error: attempted attachment {fil} does not appear to exist"
            else:
                message.attach(str(fil))
    nl = "\n"
    message.plain += f"{nl}{nl}If you have any question you can reach us at:"
    message.plain += f"{nl}clinicalgenomics@gu.se"
    message.plain += f"{nl}"
    message.plain += f"{nl}Best regards,"
    message.plain += f"{nl}Clinical Genomics Gothenburg"
    mailer.send(message)
