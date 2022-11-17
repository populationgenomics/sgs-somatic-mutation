#!/usr/bin/env python3
"""
submit job by dataproc
Suggested and programmed by Zhili<zhilizheng@outlook.com>
Tested by Zhen Qiao
License: MIT
Please read the document from analysis-runner
https://anaconda.org/cpg/analysis-runner
"""

import os
import hailtop.batch as hb
from analysis_runner import dataproc
import click

from cpg_utils.config import get_config

from cpg_utils.hail_batch import (
    authenticate_cloud_credentials_in_job,
    copy_common_env,
    remote_tmpdir,
    output_path
)

config = get_config()

service_backend = hb.ServiceBackend(
    billing_project=config['hail']['billing_project'], remote_tmpdir=remote_tmpdir()
)

batch = hb.Batch(name='Dataproc batch', backend=service_backend)

@click.command()
@click.option("--script", help="script to run")
@click.option("--jobname", help="jobname")
@click.option("--time", default='2h', help="max time")
@click.option("--nworker", default=2, help="number of worker")
@click.option("--worker", default=None, help="type of worker: e.g. 'n1-highmem-8'")
def submit(script, jobname, time, nworker, worker):
    cluster = dataproc.setup_dataproc(
            batch,
            max_age=time,
            num_workers=nworker,
            worker_machine_type=worker,
            packages=['click', 'selenium'],
            init=['gs://cpg-reference/hail_dataproc/install_common.sh'],
            cluster_name='Dataproc cluster',
            )
    cluster.add_job(script, job_name=jobname)
    # Don't wait, which avoids resubmissions if this job gets preempted.
    batch.run(wait=False)

if __name__ == "__main__":
    submit()

