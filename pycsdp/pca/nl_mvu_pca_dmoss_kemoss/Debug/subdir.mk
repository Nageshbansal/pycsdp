################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../framework.c \
../io.c \
../kernel.c \
../main.c \
../pca.c \
../report.c 

OBJS += \
./framework.o \
./io.o \
./kernel.o \
./main.o \
./pca.o \
./report.o 

C_DEPS += \
./framework.d \
./io.d \
./kernel.d \
./main.d \
./pca.d \
./report.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -fPIC -DVERSION=\"1.0\" -O0 -g3 -Wall -c -fmessage-length=0 -DNOSHORTS -ansi -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


