
import requests
url='https://www.google.com/'
response1=requests.get(url)
f = open("Google.txt",'w')
f.write(response1.text)


