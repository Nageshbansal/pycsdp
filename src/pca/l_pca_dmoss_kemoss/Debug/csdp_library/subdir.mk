################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../csdp_library/Fnorm.c \
../csdp_library/add_mat.c \
../csdp_library/addscaledmat.c \
../csdp_library/allocmat.c \
../csdp_library/calc_dobj.c \
../csdp_library/calc_pobj.c \
../csdp_library/chol.c \
../csdp_library/copy_mat.c \
../csdp_library/easysdp.c \
../csdp_library/freeprob.c \
../csdp_library/initparams.c \
../csdp_library/initsoln.c \
../csdp_library/linesearch.c \
../csdp_library/make_i.c \
../csdp_library/makefill.c \
../csdp_library/mat_mult.c \
../csdp_library/mat_multsp.c \
../csdp_library/matvec.c \
../csdp_library/norms.c \
../csdp_library/op_a.c \
../csdp_library/op_at.c \
../csdp_library/op_o.c \
../csdp_library/packed.c \
../csdp_library/psd_feas.c \
../csdp_library/qreig.c \
../csdp_library/readprob.c \
../csdp_library/readsol.c \
../csdp_library/sdp.c \
../csdp_library/solvesys.c \
../csdp_library/sortentries.c \
../csdp_library/sym_mat.c \
../csdp_library/trace_prod.c \
../csdp_library/tweakgap.c \
../csdp_library/user_exit.c \
../csdp_library/writeprob.c \
../csdp_library/writesol.c \
../csdp_library/zero_mat.c 

OBJS += \
./csdp_library/Fnorm.o \
./csdp_library/add_mat.o \
./csdp_library/addscaledmat.o \
./csdp_library/allocmat.o \
./csdp_library/calc_dobj.o \
./csdp_library/calc_pobj.o \
./csdp_library/chol.o \
./csdp_library/copy_mat.o \
./csdp_library/easysdp.o \
./csdp_library/freeprob.o \
./csdp_library/initparams.o \
./csdp_library/initsoln.o \
./csdp_library/linesearch.o \
./csdp_library/make_i.o \
./csdp_library/makefill.o \
./csdp_library/mat_mult.o \
./csdp_library/mat_multsp.o \
./csdp_library/matvec.o \
./csdp_library/norms.o \
./csdp_library/op_a.o \
./csdp_library/op_at.o \
./csdp_library/op_o.o \
./csdp_library/packed.o \
./csdp_library/psd_feas.o \
./csdp_library/qreig.o \
./csdp_library/readprob.o \
./csdp_library/readsol.o \
./csdp_library/sdp.o \
./csdp_library/solvesys.o \
./csdp_library/sortentries.o \
./csdp_library/sym_mat.o \
./csdp_library/trace_prod.o \
./csdp_library/tweakgap.o \
./csdp_library/user_exit.o \
./csdp_library/writeprob.o \
./csdp_library/writesol.o \
./csdp_library/zero_mat.o 

C_DEPS += \
./csdp_library/Fnorm.d \
./csdp_library/add_mat.d \
./csdp_library/addscaledmat.d \
./csdp_library/allocmat.d \
./csdp_library/calc_dobj.d \
./csdp_library/calc_pobj.d \
./csdp_library/chol.d \
./csdp_library/copy_mat.d \
./csdp_library/easysdp.d \
./csdp_library/freeprob.d \
./csdp_library/initparams.d \
./csdp_library/initsoln.d \
./csdp_library/linesearch.d \
./csdp_library/make_i.d \
./csdp_library/makefill.d \
./csdp_library/mat_mult.d \
./csdp_library/mat_multsp.d \
./csdp_library/matvec.d \
./csdp_library/norms.d \
./csdp_library/op_a.d \
./csdp_library/op_at.d \
./csdp_library/op_o.d \
./csdp_library/packed.d \
./csdp_library/psd_feas.d \
./csdp_library/qreig.d \
./csdp_library/readprob.d \
./csdp_library/readsol.d \
./csdp_library/sdp.d \
./csdp_library/solvesys.d \
./csdp_library/sortentries.d \
./csdp_library/sym_mat.d \
./csdp_library/trace_prod.d \
./csdp_library/tweakgap.d \
./csdp_library/user_exit.d \
./csdp_library/writeprob.d \
./csdp_library/writesol.d \
./csdp_library/zero_mat.d 


# Each subdirectory must supply rules for building sources it contributes
csdp_library/%.o: ../csdp_library/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -fPIC -DVERSION=\"1.0\" -O0 -g3 -Wall -c -fmessage-length=0 -DNOSHORTS -ansi -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


