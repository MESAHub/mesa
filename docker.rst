################
MESA with Docker
################

Docker is a tool that allows you to package an application and its dependencies in a virtual container.
This makes it easy to deploy and run applications without worrying about the underlying system configuration.

Here we provide a Docker image that contains the MESA and all its dependencies.
This image is based on the `mesa-docker <https://hub.docker.com/r/mesastar/mesa-docker>`_ repository.

Build the container
===================

To build the container, you need to have Docker installed on your system.

You can build the container by running the following command in the terminal in your MESA directory:

.. code-block:: bash

    docker build --platform linux/amd64 -t mesa-docker .

This command will build the container and tag it with the name ``mesa-docker``.

Run with Docker
===============

.. code-block:: bash

    docker run mesa-docker


Run with Singularity/Apptainer
==============================
