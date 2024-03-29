## How to put your github repo submodule in the bandframework/software directory

> First fork the bandframework repo on your local computer:
> - Go to *[bandframework repo](https://github.com/bandframework/bandframework/)* and click fork on upper right
> - Go to your github account and a copy of the bandframework repo should appear there.
> - Clone the `develop` branch of your forked repo to your local computer:

`git clone -b develop <link to bandframework forked repo>`

> To get the content of any submodules (e.g., those in `bandframework/software`) execute the following command:

`git submodule update --init --recursive`

> To add your software as a submodule

```
cd bandframework/software
git submodule add <url of your custom software repo>
  ```


> Now update your forked repo with this change
 ```
  git add .
  git commit -m "added submodule"
  git push origin develop
```

> Put your software on bandframework develop branch:
>- Go to the *[develop branch of the bandframework](https://github.com/bandframework/bandframework/tree/develop)*
>- Click on the *[pull request tab](https://github.com/bandframework/bandframework/pulls)* at the very top.
>- Click on the green `New pull request` button.
>- Click compare across all.
>- Follow the instructions and finally submit the pull request with a comment for the approval. 

### Later when you have made changes to your software and want to update the bandframework with this change do the following

> While on the develop branch (or a branch created from develop) go into bandframework/software and issue
`  git submodule update --remote`
> Then, add the changes for your_submodule: 
 ```
 git add your_submodule
 git commit -m "updating the submodule your_submodule"
 git push origin develop
 ```

## BAND-recommended way to log data for BAND packages added via submodules

> Software directly hosted in this repository has data logging enabled to track number of clones, etc.
> for a longer time period than github's.
>
> To track this data for packages hosted elsewhere, we recommend using https://github.com/marketplace/actions/github-repo-stats
>
>  Users should follow the tutorial https://github.com/jgehrcke/github-repo-stats/wiki/Tutorial
