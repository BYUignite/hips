docker build -t hips_img .

docker run --rm --name hips -it -v $(pwd)/..:/app/hips hips_img   