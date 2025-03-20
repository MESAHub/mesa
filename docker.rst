################
MESA with Docker
################

Docker is a tool that allows you to package an application and its dependencies in a virtual container.
This makes it easy to deploy and run applications without worrying about the underlying system configuration.

We provide a Docker image that contains the MESA and all its dependencies.


Obtain the Docker Image
=======================

You can obtain the Docker image from the Docker Hub repository by running the following command:

.. code-block:: bash

    docker pull mesastar/mesa-docker:XX.XX.X


Run with Docker
===============

You can run the Docker container by running the following command:

.. code-block:: bash

    docker run mesa-docker


Run with Singularity/Apptainer
==============================


Build the Container (for Developers)
====================================

To build the container, you need to have Docker installed on your system.

You can build the container by running the following command in the terminal in your MESA directory:

.. code-block:: bash

    docker build --platform linux/arm64 -t mesa-docker .

This command will build the container and tag it with the name ``mesa-docker``.
