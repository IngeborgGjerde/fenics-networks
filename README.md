# Implementing network models in FEniCS

This is a repo for experimenting with using `networkx` to make graph meshes and solve equations on these in `FEniCS`. 

## Installation
We use `networkx` combined with the mixed-dimensional library [`fenics_ii`](https://github.com/MiroK/fenics_ii) created by Miroslav Kuchta. 

The environment is provided as a docker container. The container can be built and run locally by executing

```shell
git clone https://github.com/IngeborgGjerde/fenics-networks/
cd fenics-networks
docker build --no-cache -t fenics-networks .
docker run --name networkfenics -v "$(pwd):/home/fenics/shared" -d -p 127.0.0.1:8888:8888 fenics-networks 'jupyter-notebook --ip=0.0.0.0'
```

You can then enter the container by running 
```shell
docker exec -it networkfenics /bin/bash
```
To connect it to a jupyter notebook, run
```shell
docker logs networkfenics
```
and enter the html-links it provides in your browser.

[<img alt="alt_text" width="400px" src="pial-network.png" />]([https://www.google.com/](https://github.com/IngeborgGjerde/fenics-networks/pial-network.png?raw=true))


## Citing

The manuscript corresponding to this repo will appear in the not-so-distant future. Feel free to contact me if you would like to read the manuscript. If you use this repo for your own work you can cite the placeholder *Network model for perivascular fluid flow driven by vasomotion* (in preparation) by Ingeborg Gjerde, Marie Rognes and Barbara Wohlmuth. 
