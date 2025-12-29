#!/usr/bin/env python3

import glob

import matplotlib.pyplot as plt
import mesa_reader as mr
from matplotlib.animation import FFMpegWriter, FuncAnimation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from static_HISTORY_check import (MesaView, read_header_columns,
                                  setup_hr_diagram_params)


def make_cmd_rotation_video(
    history_file="../LOGS/history.data",
    outfile="cmd_interp_rad_rotation.mp4",
    fps=30,
    total_seconds=10,
):
    md = mr.MesaData(history_file)
    md = MesaView(md, 5)

    all_cols, filter_columns = read_header_columns(history_file)
    hr_color, hr_mag, hr_xlabel, hr_ylabel, color_index = setup_hr_diagram_params(
        md, filter_columns
    )

    interp_rad = md.Interp_rad

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    sc = ax.scatter(hr_color, hr_mag, interp_rad, s=5)

    ax.set_xlabel(hr_xlabel)
    ax.set_ylabel(hr_ylabel)
    ax.set_zlabel("Interp_rad")
    ax.invert_yaxis()

    nframes = fps * total_seconds
    hold_frames = fps * 1  # 1 second showing plain CMD

    # Start looking straight down the Interp_rad axis: CMD view
    start_elev = 90
    start_azim = -90

    # End at an oblique 3D view
    end_elev = 5
    end_azim = 270

    def init():
        ax.view_init(elev=start_elev, azim=start_azim)
        return (sc,)

    def update(frame):
        # hold at start
        if frame <= hold_frames:
            elev = start_elev
            azim = start_azim

        # hold at end
        elif frame >= nframes - hold_frames:
            elev = end_elev
            azim = end_azim

        else:
            # transition region
            t = (frame - hold_frames) / (nframes - 2 * hold_frames)
            elev = start_elev + t * (end_elev - start_elev)
            azim = start_azim + t * (end_azim - start_azim)

        ax.view_init(elev=elev, azim=azim)
        return (sc,)

    anim = FuncAnimation(
        fig,
        update,
        init_func=init,
        frames=nframes,
        interval=1000 / fps,
        blit=False,
    )

    writer = FFMpegWriter(fps=fps, bitrate=2000)
    anim.save(outfile, writer=writer)


def main():
    history_file = glob.glob("../LOGS/history.data")[0]
    make_cmd_rotation_video(history_file)


if __name__ == "__main__":
    main()
