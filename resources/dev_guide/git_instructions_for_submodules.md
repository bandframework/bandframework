## How to put your github repo submodule in the bandframework/software directory

> First fork the bandframework repo on your local computer:
> - Go to *[bandframework repo](https://github.com/bandframework/bandframework/)* and click fork on upper right
> - Go to your github account and a copy of the bandframework repo should appear there.
> - Clone the `development` branch of your forked repo to your local computer:

`git clone -b development <link to bandframework forked repo>`

> To get the content of any submodules (e.g., those in `bandframework/software`) execute the following command:

`git submodule update --init --recursive`

> To add your software as a submodule

```
cd bandframework/software
git submodule add <url of your custom softwear repo>
  ```


> Now update your forked repo with this change
 ```
  git add .
  git commit -m "added submodule"
  git push origin development
```

> Put your software on bandframework development branch:
>- Go to the *[development branch of the bandframework](https://github.com/bandframework/bandframework/tree/development)*
>- Click on the *[pull request tab](https://github.com/bandframework/bandframework/pulls)* at the very top.
>- Click on the green `New pull request` button.
>- Click compare accross all.
>- Follow the instructions and finally submit the pull request with a comment for the approval. 

### Later when you have made changes to your software and want to update the bandframework with this change do the following

> Go into the bandframework/software/your_sub_module
`  git pull origin`
> Go back to the previous directory: bandframework/software
 ```
 git add .
 git commit -m "updating the submodule"
 git push origin development
 ```
 > Then you have to repeat the above procedure under "Put your software apear on bandframework development branch"
