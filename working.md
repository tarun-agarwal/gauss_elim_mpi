**Ryan Ordille** (260399372)

ryan.ordille@mail.mcgill.ca

ECSE 420: Parallel Computing, Fall 2014

Lab 2 - 07 November, 2014

# Setup

Apple MacBook Air, 13", Mid 2011 (model 4,2)

**Operating System**: OS X 10.10 Yosemite

**Processor**: 1.7 GHz Intel Core i5, 2 cores, 256 KB L2 cache per core, 3 MB L3 cache

**Other Information**: 4 GB 1333 MHz DDR3 RAM, Intel HD Graphics 3000 384 MB, SSD drive

Using `open-mpi-1.8.3`

# 1. Gaussian Elimination

## 1.1 Parallelization Strategy

## 1.2 Serial Implementation

See `GESerial.c`.

### Results

Serial implementation output:

![](./pictures/serial.png)

Serial implementation timing:

![](./pictures/serial_time.png)

    [Serial] (1024x1024) 1.834153e+00 s
    [Serial] (2048x2048) 1.491354e+01 s
    [Serial] (4096x4096) 1.204770e+02 s

# 2. Message Passing Program

## 2.1 Point-to-Point

### 2.1.1 Continuous

### 2.1.1 Circular

### Results

## 2.2 Broadcast/Collective

### 2.2.1 Continuous

### 2.2.2. Circular

### Results

# 3. Experiment


