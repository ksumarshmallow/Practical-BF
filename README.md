 **🧬Practical-BF🧬**
================

## **Описание**
В этом репозитории собраны решения домашних заданий по дисциплине "Практическая биоинформатика", которые были выполнены в ходе обучения.

-------------

##  **Примерная структура проекта**
```markdown
Practical-BF/
├── HW1/
│   ├── task1/
│   │   └── task1.py
│   ├── task2/
│   │   └── task2.py
│   ├── HW1_Malkova.ipynb
│   └── HW1_Malkova.html
├── HW2/
│   ├── task1/
│   │   └── task1.py
│   ├── task2/
│   │   └── task2.py
│   ├── HW2_Malkova.ipynb
│   └── HW2_Malkova.html
├── ...
├── utils/
│   ├── run_cmd.py
│   └── ...
├── README.md
├── pyproject.toml
└── environment.yaml
```

* В каждой папке с домашним заданием (**`HW1`**, **`HW2`**, ...) содержит папки с названиями заданий (**`task1`**, **`task2`**, ...)
* Каждая папка с заданиями содержит `.py`-файл(-ы), которые были использованы для решения соответствующего задания
* Также каждая папка с домашним заданием содержит `.ipynb`-файл, который совмещает решения всех заданий в данной домашней работе. Также в нем приведены идея и выводы по каждому из заданий (если это было необходимо)
* Для каждого `.ipynb`-файла есть соответствующий `.html`-файл. Для его просмотра можно:
    1. Скачать `.html` и открыть локально. Но на телефоне могут быть проблемы
    2. К ссылке с `.html`, лежащей в гит-репозитории, добавить приставку **https://html-preview.github.io/?url=**. Например, [вот ссылка](https://html-preview.github.io/?url=https://github.com/ksumarshmallow/Practical-BF/blob/main/HW1/HW1_Malkova.html) на `.html` первого ДЗ.
* Директория **`utils`** содержит файлы с различными вспомогательными функциями, которые могут быть полезны при решении заданий (вне зависимости от номера домашней работы)

------------
## **Установка**

### 1. Клонирование репозитория
```bash
cd ~
git clone https://github.com/ksumarshmallow/Practical-BF.git
cd Practical-BF
```

### 2. Создание виртуальной среды

Можете использовать любую виртуальную среду. Здесь представлены два варианта по их установке: виртуальная среда Python и conda

#### 2.1 Виртуальная среда Python

0. Установите необходимую версию python (можете попробовать использовать версию по умолчанию. Но в этом проекте рекомендуется использовать python 3.12). Для установки python определенной версии без прав **sudo** можно использовать `pyenv`:

```bash
cd ~
curl https://pyenv.run | bash
source ~/.bashrc
pyenv install 3.12.4
```

1. Создайте виртуальную среду:

    ```bash
    python3.12 -m venv .venv
    ```

2. Активируйте виртуальную среду:

    - Для Linux/MacOS:

      ```bash
      source .venv/bin/activate
      ```

    - Для Windows:

      ```bash
      .venv\Scripts\activate
      ```

#### 2.2 Виртуальная среда conda

1. Создайте виртуальную среду conda с желаемой версией python (здесь использовался 3.12.4)
    ```bash
    conda create -n practical_bf python=3.12.4
    ```

    Можете создать виртуальную среду с вресией python по умолчанию:
   ```bash
    conda create -n practical_bf
    ``` 

3. Активируйте виртуальную среду:
    ```bash
    conda activate practical_bf
    ```

4. Установите в виртаульной среде ipykernel - необходимо для работы в .ipynb в той же среде
    ```bash
    conda install ipykernel
    ipython kernel install --user --name=practical_bf
    ```

### 3. Установка зависимостей
```bash
pip install -r requirements.txt
```

### 4. Установка инструментов

#### Если у вас есть права `sudo`
1. **Установка зависимостей**
```bash
sudo apt-get update
sudo apt-get install build-essential autoconf zlib1g-dev
```

2. **Установка инструментов** (только последняя версия)
```bash
sudo apt-get install bedtools samtools bcftools
```

#### Установка `samtools`

1. **Скачайте исходный код с официального сайта**. Выберите нужную версию по [этой ссылке](https://github.com/samtools/samtools/releases/) и замените версию в команде:

    ```bash
    wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2
    tar -xjf samtools-1.20.tar.bz2
    cd samtools-1.20
    ```

2. **Создайте пользовательский каталог для установки**:

    ```bash
    mkdir -p $HOME/local
    ```

3. **Соберите и установите `samtools` в пользовательский каталог**:

    ```bash
    ./configure --prefix=$HOME/local
    make
    make install
    ```

4. **Добавьте каталог в `PATH`**:

    ```bash
    export PATH=$HOME/local/bin:$PATH
    source ~/.bashrc
    ```

5. **Проверьте установку**:

    ```bash
    samtools --version
    ```

#### Установка `bedtools`

1. **Скачайте исходный код с официального сайта**. Выберите нужную версию по [ссылке](https://bedtools.readthedocs.io/en/latest/content/history.html) и замените версию в команде:

    ```bash
    wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz
    tar -xzf bedtools-2.31.1.tar.gz
    cd bedtools2
    ```

2. **Соберите `bedtools`**:

    ```bash
    make
    ```

3. **Создайте каталог для установки и скопируйте исполняемый файл `bedtools` в него**:

    ```bash
    mkdir -p $HOME/.local/bin
    cp bin/bedtools $HOME/.local/bin/
    ```

4. **Добавьте каталог в `PATH`**:

    ```bash
    export PATH="$HOME/.local/bin:$PATH"
    source ~/.bashrc
    ```

5. **Проверьте установку**:

    ```bash
    bedtools --version
    ```

#### Установка `bcftools`

1. **Скачайте исходный код с официального сайта**. Выберите нужную версию по [ссылке](https://github.com/samtools/bcftools/releases/) и замените версию в команде:

    ```bash
    wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2
    tar -xjf bcftools-1.19.tar.bz2
    cd bcftools-1.19
    ```

2. **Соберите и установите `bcftools` в пользовательский каталог**:

    ```bash
    make
    mkdir -p $HOME/.local/bin
    cp bcftools $HOME/.local/bin/
    ```

3. **Добавьте каталог в `PATH`**:

    ```bash
    export PATH="$HOME/.local/bin:$PATH"
    source ~/.bashrc
    ```

4. **Проверьте установку**:

    ```bash
    bcftools --version
    ```

------------
# **Авторы**
*Малкова Ксения Эдуардовна*, студент группы (я еще не знаю)