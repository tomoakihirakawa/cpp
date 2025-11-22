#!/bin/bash

# リモートリポジトリから最新の変更を取得
git fetch origin

# mainブランチに切り替え、最新状態に更新
git checkout main
git pull origin main

# testブランチをmainにマージ
git merge test

# mainブランチの変更をリモートにプッシュ
git push origin main

# 再びtestブランチに切り替え
git checkout test
