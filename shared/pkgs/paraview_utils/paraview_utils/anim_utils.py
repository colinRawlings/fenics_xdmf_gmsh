"""
Some helpful code for exporting animations from paraview

In particular hopefully avoids the frequent crashes with more
calls to Render
"""

#############################################################################
# Imports
#############################################################################

import subprocess as sp

# add python paraview bin to python path so that the import works(?) probably a terrible idea
# as would see all paraview packages ...

from paraview import simple as pv

#############################################################################
# Defintions
#############################################################################

FRAME_FILE_EXTENSION = ".jpg"
MOVIE_FILE_EXTENSION = ".mov"

###############################################################
# Functions
###############################################################


def create_animation(output_prefix="output",
                     init_camera_position=(1, 1, 1),
                     DAzimuth=90,
                     DElevation=20,
                     DDolly=1.5):

    # init

    sp.check_output(['where', 'avconv'])

    # set abs az = 45, abs el = 45

    camera = pv.GetActiveCamera()
    camera.SetPosition(*init_camera_position)
    camera.SetFocalPoint(0, 0, 0)
    camera.SetViewUp(0, 0, 1)
    pv.Render()

    # n steps

    anim = pv.GetAnimationScene()
    num_steps = int(anim.EndTime) + 1

    del_az = DAzimuth / num_steps
    del_elv = DElevation / num_steps
    del_dolly = DDolly**(1 / num_steps)

    # create animation

    anim.GoToFirst()
    pv.Render()

    for time_index in range(num_steps):

        print("Processing frame: {} of {}".format(time_index + 1, num_steps))
        
        pv.SaveScreenshot(output_prefix + "_frame{}{}".format(time_index, FRAME_FILE_EXTENSION))
        pv.Render()

        anim.GoToNext()
        pv.Render()

        # update camera

        camera.Azimuth(del_az)
        pv.Render()

        camera.Elevation(del_elv)
        pv.Render()

        camera.Dolly(del_dolly)
        pv.Render()

        #

    # video

    cmd_list = [
        "avconv", "-i", output_prefix + "_frame%d" + FRAME_FILE_EXTENSION, "-qscale", "1",
        output_prefix + MOVIE_FILE_EXTENSION
    ]

    proc = sp.Popen(cmd_list, stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = proc.communicate()

    assert proc.returncode == 0, "Movie creation failed.{}{}".format(
        stdout, stderr)
