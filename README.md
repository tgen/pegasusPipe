# Notes while installing on Dback

My current thoughts for organizing this project better:

So far, I've refactored the filesystem links in centralPipe and pegasusPipe to be relative to the install location. For example, Pegasus looks for job scripts in "<install_dir>/jobScripts/" instead of "/home/tgenjetstream/pegasusPipe/jobScripts/". This means that we can install all the modules (centralPipe, pegasusPipe, medusaPipe) as a single bundle, Here is how I've installed it in my home, I'm calling it jetstream:

```
[rrichholt@dback-login1]~$ ls jetstream/
centralPipe/  constants/  pegasusPipe/  tests/
```

Since constants are used by all the other modules, they should really be considered a separate module. So, I've moved the constants directory out of central pipe and under the main jetstream/ directory.

# Big differences between slurm and torque

First job that I tried to launch fails without any errors, and no oeFiles. There is a problem with the way Slurm is handling the --output and --error directives:

To test this I've made a small job script and set out to be ~/wtf.output: 

```shell
[rrichholt@dback-login1]~/scripts$ cat oeTest.sh
#!/usr/bin/env bash
#SBATCH --output ~/wtf.output
echo hello world, this is ${name}
```

Now I submit the script:

```shell
[rrichholt@dback-login1]~/scripts$ sbatch oeTest.sh --export ALL,name=Jeff
Submitted batch job 3978
```

And after a couple minutes there is still no output file. If I look up the job with `scontrol`, you can see that StdOut is set to an invalid location. 

```shell
rrichholt@dback-login1]~/scripts$ scontrol show job 3978
JobId=3978 JobName=oeTest.sh
   UserId=rrichholt(735597512) GroupId=Domain Users(1689184361) MCS_label=N/A
   Priority=888397 Nice=0 Account=tgen-204000 QOS=normal
   JobState=FAILED Reason=NonZeroExitCode Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=1:0
   RunTime=00:00:01 TimeLimit=14-00:00:00 TimeMin=N/A
   SubmitTime=2017-09-11T09:50:39 EligibleTime=2017-09-11T09:50:39
   StartTime=2017-09-11T09:50:39 EndTime=2017-09-11T09:50:40 Deadline=N/A
   PreemptTime=None SuspendTime=None SecsPreSuspend=0
   Partition=defq AllocNode:Sid=dback-login1:5290
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=dback-c1-n01
   BatchHost=dback-c1-n01
   NumNodes=1 NumCPUs=1 NumTasks=0 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=1,mem=4500M,node=1
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryCPU=4500M MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   Gres=(null) Reservation=(null)
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=/home/rrichholt/scripts/oeTest.sh --export ALL,name=Jeff
   WorkDir=/home/rrichholt/scripts
   StdErr=/home/rrichholt/scripts/~/wtf.output
   StdIn=/dev/null
   StdOut=/home/rrichholt/scripts/~/wtf.output
   Power=
```

Okay, so maybe glob expansion isn't applied? Now I'll try with an absolute path as the --output location: /home/rrichholt/wtf.output...

No console logs here, same as above, I just changed the SBATCH directive.

And it works as expected with `scontrol` showing `StdOut=/home/rrichholt/wtf.output`. 

Now let's try with a variable name:

```shell
[rrichholt@dback-login1]~/scripts$ cat oeTest.sh
#!/usr/bin/env bash
#SBATCH --output /home/rrichholt/${name}.output
echo hello world, this is ${name}
```

And now this is weird:

```shell
[rrichholt@dback-login1]~/scripts$ scontrol show job 3981
JobId=3981 JobName=oeTest.sh
   UserId=rrichholt(735597512) GroupId=Domain Users(1689184361) MCS_label=N/A
   Priority=888929 Nice=0 Account=tgen-204000 QOS=normal
   JobState=COMPLETED Reason=None Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:00:00 TimeLimit=14-00:00:00 TimeMin=N/A
   SubmitTime=2017-09-11T10:00:39 EligibleTime=2017-09-11T10:00:39
   StartTime=2017-09-11T10:00:39 EndTime=2017-09-11T10:00:39 Deadline=N/A
   PreemptTime=None SuspendTime=None SecsPreSuspend=0
   Partition=defq AllocNode:Sid=dback-login1:5290
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=dback-c1-n01
   BatchHost=dback-c1-n01
   NumNodes=1 NumCPUs=1 NumTasks=0 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=1,mem=4500M,node=1
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryCPU=4500M MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   Gres=(null) Reservation=(null)
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=/home/rrichholt/scripts/oeTest.sh
   WorkDir=/home/rrichholt/scripts
   StdErr=/home/rrichholt/${name}.output
   StdIn=/dev/null
   StdOut=/home/rrichholt/${name}.output
   Power=
```

The variable never gets expanded, it writes to that literal location:

```shell
[rrichholt@dback-login1]~/scripts$ ls ~/
bin/  jetstream/  ${name}.output  projects/  scripts/
```
