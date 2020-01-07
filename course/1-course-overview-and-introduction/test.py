#!/usr/bin/env python3

import datetime
import gantt

gantt.define_font_attributes(fill='black',
                             stroke='black',
                             stroke_width=0,
                             font_family="Arial")

gantt.add_vacations(datetime.date(2020,1,20))
gantt.add_vacations(datetime.date(2020,2,17))
gantt.add_vacations(datetime.date(2020,3,9),datetime.date(2020,3,14))
gantt.add_vacations(datetime.date(2020,4,23));


# Create a task
tintro = gantt.Task(name='Course Intro',
                    start=datetime.date(2020, 1, 7), duration=1)
tproblem = gantt.Task(name='Problem Statement',
                      start=datetime.date(2020, 1, 8), duration=5)
tcv = gantt.Task(name='Computer vision',
                 start=datetime.date(2020,1,13),duration=5,
                 depends_of=[tproblem])
tnerf = gantt.Task(name='Nerf gun analysis',
                   start=datetime.date(2020,1,20),duration=5)
tfiring = gantt.Task(name='Firing circuit development',
                     start=datetime.date(2020,1,30),duration=10,
                     depends_of=[tnerf])
tturret = gantt.Task(name='Turret subsystem',
                     start=datetime.date(2020,2,10),duration=10)
tmodeling = gantt.Task(name='Turret system modeling',
                       start=datetime.date(2020,3,2),duration=5,
                       depends_of=[tturret])
tcontrol = gantt.Task(name='Turret position control',
                      start=datetime.date(2020,3,16),duration=10,
                      depends_of=[tmodeling])
tballistics = gantt.Task(name='Ballistics',
                         start=datetime.date(2020,3,25),duration=10,
                         depends_of=[tcontrol,tfiring])
tintegration = gantt.Task(name='System integration',
                          start=datetime.date(2020,4,9),duration=9,
                          depends_of=[tcv,tfiring,tcontrol,tballistics])

mfiring = gantt.Milestone(name='Firing subsystem demo',start=datetime.date(2020,4,13),depends_of=[tfiring])
mfinal = gantt.Milestone(name='Final demo',start=datetime.date(2020,4,27),
                         depends_of=[tintegration])

# Create a project
pEW309 = gantt.Project(name='Sp2020 EW309 0111 0311');
pEW309.add_task(tintro)
pEW309.add_task(tproblem)
pEW309.add_task(tcv)
pEW309.add_task(tnerf)
pEW309.add_task(tfiring)
pEW309.add_task(tturret)
pEW309.add_task(tmodeling)
pEW309.add_task(tcontrol)
pEW309.add_task(tballistics)
pEW309.add_task(tintegration)
pEW309.add_task(mfiring)
pEW309.add_task(mfinal)

# render
pEW309.make_svg_for_tasks(filename='sp2020-ew309-gantt.svg',
                          today=datetime.date(2020,1,7))
