qsub \
    -t 1-$(cat /path1/to/samples.joint_genotyping.list | wc -l) \
    -q short.qc \
    /path2/to/submit_jobs.joint_genotyping.sh /path1/to/samples.joint_genotyping.list
