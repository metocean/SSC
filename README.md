# SSC
Sustainable Sea project

Build the docker
====

Build the image from the dockerfile:

```
$ cd /where/you/save/Dockerfile
$ sudo docker build -t metocean/ssc .
```

After this you will have an image called `ssc` which contains:
- the schism model compile ( you can type schism from anywhere inside)
- the combine_output10 program to compile the output
- the git hub folder with all the code
everything will be in /home/user/SSC

Create a container
====

Once you have set up the image you need to run it ( that will create a container)
```
$ docker run --name SSC -it metocean/ssc
```
or with mounting a existing folder
```
$ docker run --name SSC -it -v /home/ross/mycode:/home/user/ metocean/ssc
```

Run the Wrapper
====

Once inside you can go in /home/user/SSC and do:
```
$ git pull to download the latest code
```
To run my code you can test with
```
$ python /home/user/SSC/func/SSC.py --yaml /home/user/SSC/test.yaml
```
or insert it into another function:
```
from SSC import wrapper
wrapper(parameters)
```

Save the work to the cloud
====

add the library you have installed in the Dockerfile
update the image from the container:
```
$ docker commit -m "new ssc_tide version" [container id] metocean/ssc_tide:latest
```
where you can find the container id with:
```
$ sudo docker ps -a
```


Save the work locally
====

if you want to export smth out of the container you can do
````
$ sudo docker cp SSC:/home/user/result/file1.nc /where/you/want/
````

Exiting the container and login back in
====

if you exit the container you can go back inside with
````
$ sudo docker start SSC
$ sudo docker attach SSC
````
