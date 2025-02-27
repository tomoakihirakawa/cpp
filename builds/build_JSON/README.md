# Contents
- [🐋 JSONクラス](#-jsonクラス)
    - [⛵ ⛵ C++でのJSON操作に関する実装と使用方法](#--cでのjson操作に関する実装と使用方法)
        - [🪼 🪼 前提条件](#--前提条件)
        - [🪼 🪼 JSONクラス](#--jsonクラス)
            - [🪸 🪸 コンストラクタ](#--コンストラクタ)
            - [🪸 🪸 操作](#--操作)
        - [🪼 🪼 使用例](#--使用例)
        - [🪼 🪼 ファイル出力](#--ファイル出力)
        - [🪼 🪼 その他の機能](#--その他の機能)
        - [🪼 🪼 注意点](#--注意点)
        - [🪼 🪼 todo](#--todo)


---
# 🐋 JSONクラス 

## ⛵ ⛵ C++でのJSON操作に関する実装と使用方法  

### 🪼 🪼 前提条件  

このドキュメントでは、`basic.hpp` と呼ばれる基本的なヘッダーファイルが含まれていると仮定します。また、サンプルのJSONデータが `./sample.json` に格納されていると仮定します。

### 🪼 🪼 JSONクラス  

この実装には、`JSON`という名前のC++クラスがあります。このクラスは、内部的に`std::map<std::string, std::vector<std::string>>`を持っており、JSONオブジェクトのデータを保存します。

#### 🪸 🪸 コンストラクタ  

このクラスにはいくつかのコンストラクタがあります:

1. ファイル名を引数として取る
```cpp
JSON json("./sample.json");
```
2. `std::ifstream` オブジェクトを引数として取る
```cpp
JSON json(std::ifstream("./sample.json"));
```

#### 🪸 🪸 操作  

- キーを用いた値の取得
```cpp
std::cout << json["translate"] << std::endl;
```
- キーを用いた値の設定
```cpp
json["price"] = {"10."};
```

### 🪼 🪼 使用例  

以下は、このクラスの簡単な使用例です。

```cpp
{
std::cout << magenta << "1. ファイル名でJSONをコンストラクト" << colorReset << std::endl;
JSON json("./sample.json");
std::cout << json["translate"] << std::endl;
}
```

```cpp
{
std::cout << red << "2. ifstreamでJSONをコンストラクト" << colorReset << std::endl;
JSON json(std::ifstream("./sample.json"));
json["price"] = {"10."};
}
```

### 🪼 🪼 ファイル出力  

JSONオブジェクトをファイルに出力するには、`operator<<`を使用します。

```cpp
std::ofstream os("./output.json");
os << json;
os.close();
```

### 🪼 🪼 その他の機能  

- `find` メソッドを使い、キーが存在するかどうかを調べることができます。
- `contains` メソッドも同様の機能を提供します。

### 🪼 🪼 注意点  

1. キーが存在しない場合は、例外がスローされます。
2. 数値や真偽値は、`stod`や`stob`のような関数を使って変換する必要があります。

このドキュメントはC++でのJSON操作の基本的な概要を簡単に説明したものです。実際の使用ケースや要件に応じて、コードは適宜調整してください。

### 🪼 🪼 todo  

📝 以下のようなJSONデータを読み込めるようにする.

```
"mooring": ["simple_mooring",[10, 0.1, 0],[10, 0.1, 0]]
```

📝 オブジェクトリテラルに対応する．

```
"name": { "asda": 1 }
```
[../../include/basic.hpp#L2046](../../include/basic.hpp#L2046)

[./main.cpp#L6](./main.cpp#L6)

---
