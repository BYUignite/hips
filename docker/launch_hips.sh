docker build -t hips_img .   

docker run -it -v $(pwd)/..:/app/hips hips_img   