//------------------------------------------------------------------------------
// <auto-generated>                                                             
//     This code was generated by a tool.                                       
//     Changes to this file may cause incorrect behavior and will be lost if    
//     the code is regenerated.                                                 
// </auto-generated>                                                            
//------------------------------------------------------------------------------
#pragma warning disable 436
#pragma warning disable 162
#pragma warning disable 1591
using System;
using Microsoft.Quantum.Core;
using Microsoft.Quantum.Intrinsic;
using Microsoft.Quantum.Simulation.Core;

[assembly: CallableDeclaration("{\"Kind\":{\"Case\":\"Operation\"},\"QualifiedName\":{\"Namespace\":\"Estimate\",\"Name\":\"TestUnifPrime\"},\"Attributes\":[{\"TypeId\":{\"Case\":\"Value\",\"Fields\":[{\"Namespace\":\"Microsoft.Quantum.Core\",\"Name\":\"EntryPoint\",\"Range\":{\"Case\":\"Value\",\"Fields\":[{\"Item1\":{\"Line\":1,\"Column\":2},\"Item2\":{\"Line\":1,\"Column\":12}}]}}]},\"TypeIdRange\":{\"Case\":\"Value\",\"Fields\":[{\"Item1\":{\"Line\":1,\"Column\":2},\"Item2\":{\"Line\":1,\"Column\":12}}]},\"Argument\":{\"Item1\":{\"Case\":\"UnitValue\"},\"Item2\":[],\"Item3\":{\"Case\":\"UnitType\"},\"Item4\":{\"IsMutable\":false,\"HasLocalQuantumDependency\":false},\"Item5\":{\"Case\":\"Value\",\"Fields\":[{\"Item1\":{\"Line\":1,\"Column\":12},\"Item2\":{\"Line\":1,\"Column\":14}}]}},\"Offset\":{\"Item1\":38,\"Item2\":4},\"Comments\":{\"OpeningComments\":[],\"ClosingComments\":[]}}],\"Modifiers\":{\"Access\":{\"Case\":\"DefaultAccess\"}},\"SourceFile\":\"C:\\\\Users\\\\vivie\\\\github\\\\QW_Burgers\\\\Estimate\\\\Program.qs\",\"Position\":{\"Item1\":39,\"Item2\":4},\"SymbolRange\":{\"Item1\":{\"Line\":1,\"Column\":11},\"Item2\":{\"Line\":1,\"Column\":24}},\"ArgumentTuple\":{\"Case\":\"QsTuple\",\"Fields\":[[]]},\"Signature\":{\"TypeParameters\":[],\"ArgumentType\":{\"Case\":\"UnitType\"},\"ReturnType\":{\"Case\":\"UnitType\"},\"Information\":{\"Characteristics\":{\"Case\":\"EmptySet\"},\"InferredInformation\":{\"IsSelfAdjoint\":false,\"IsIntrinsic\":false}}},\"Documentation\":[]}")]
[assembly: SpecializationDeclaration("{\"Kind\":{\"Case\":\"QsBody\"},\"TypeArguments\":{\"Case\":\"Null\"},\"Information\":{\"Characteristics\":{\"Case\":\"EmptySet\"},\"InferredInformation\":{\"IsSelfAdjoint\":false,\"IsIntrinsic\":false}},\"Parent\":{\"Namespace\":\"Estimate\",\"Name\":\"TestUnifPrime\"},\"Attributes\":[],\"SourceFile\":\"C:\\\\Users\\\\vivie\\\\github\\\\QW_Burgers\\\\Estimate\\\\Program.qs\",\"Position\":{\"Item1\":39,\"Item2\":4},\"HeaderRange\":{\"Item1\":{\"Line\":1,\"Column\":11},\"Item2\":{\"Line\":1,\"Column\":24}},\"Documentation\":[]}")]
#line hidden
namespace Estimate
{
    [SourceLocation("C:\\Users\\vivie\\github\\QW_Burgers\\Estimate\\Program.qs", OperationFunctor.Body, 40, -1)]
    public partial class TestUnifPrime : Operation<QVoid, QVoid>, ICallable
    {
        public TestUnifPrime(IOperationFactory m) : base(m)
        {
        }

        String ICallable.Name => "TestUnifPrime";
        String ICallable.FullName => "Estimate.TestUnifPrime";
        public static EntryPointInfo<QVoid, QVoid> Info => new EntryPointInfo<QVoid, QVoid>(typeof(TestUnifPrime));
        protected Allocate Allocate__
        {
            get;
            set;
        }

        protected Release Release__
        {
            get;
            set;
        }

        protected ICallable Microsoft__Quantum__Arrays__Zipped
        {
            get;
            set;
        }

        protected IUnitary<Qubit> Microsoft__Quantum__Intrinsic__X
        {
            get;
            set;
        }

        protected IUnitary<(IQArray<Qubit>,IQArray<Qubit>)> AmplitudeTransduction__UnifPrime
        {
            get;
            set;
        }

        protected ICallable<IQArray<Qubit>, QVoid> ResetAll__
        {
            get;
            set;
        }

        protected ICallable<Qubit, QVoid> Reset__
        {
            get;
            set;
        }

        public override Func<QVoid, QVoid> __Body__ => (__in__) =>
        {
#line 42 "C:\\Users\\vivie\\github\\QW_Burgers\\Estimate\\Program.qs"
            var nPrecisionBits = 3L;
#line hidden
            {
#line 43 "C:\\Users\\vivie\\github\\QW_Burgers\\Estimate\\Program.qs"
                var dataRegister = Allocate__.Apply(nPrecisionBits);
#line hidden
                bool __arg1__ = true;
                try
                {
#line hidden
                    {
#line 44 "C:\\Users\\vivie\\github\\QW_Burgers\\Estimate\\Program.qs"
                        var referenceRegister = Allocate__.Apply(nPrecisionBits);
#line hidden
                        bool __arg2__ = true;
                        try
                        {
#line hidden
                            {
#line 45 "C:\\Users\\vivie\\github\\QW_Burgers\\Estimate\\Program.qs"
                                var flagQubit = Allocate__.Apply();
#line hidden
                                bool __arg3__ = true;
                                try
                                {
#line 48 "C:\\Users\\vivie\\github\\QW_Burgers\\Estimate\\Program.qs"
                                    var classicalValues = (IQArray<Boolean>)new QArray<Boolean>(false, false, true);
#line 49 "C:\\Users\\vivie\\github\\QW_Burgers\\Estimate\\Program.qs"
                                    foreach (var (classicalValue,dataQubit) in Microsoft__Quantum__Arrays__Zipped.Apply<IQArray<(Boolean,Qubit)>>((classicalValues, dataRegister)))
#line hidden
                                    {
#line 50 "C:\\Users\\vivie\\github\\QW_Burgers\\Estimate\\Program.qs"
                                        if (classicalValue)
                                        {
#line 51 "C:\\Users\\vivie\\github\\QW_Burgers\\Estimate\\Program.qs"
                                            Microsoft__Quantum__Intrinsic__X.Apply(dataQubit);
                                        }
                                    }

#line 56 "C:\\Users\\vivie\\github\\QW_Burgers\\Estimate\\Program.qs"
                                    AmplitudeTransduction__UnifPrime.Apply((QArray<Qubit>.Add(new QArray<Qubit>(flagQubit), referenceRegister), dataRegister));
#line 73 "C:\\Users\\vivie\\github\\QW_Burgers\\Estimate\\Program.qs"
                                    ResetAll__.Apply(dataRegister);
#line 73 "C:\\Users\\vivie\\github\\QW_Burgers\\Estimate\\Program.qs"
                                    ResetAll__.Apply(referenceRegister);
#line 73 "C:\\Users\\vivie\\github\\QW_Burgers\\Estimate\\Program.qs"
                                    Reset__.Apply(flagQubit);
                                }
#line hidden
                                catch
                                {
                                    __arg3__ = false;
                                    throw;
                                }
#line hidden
                                finally
                                {
                                    if (__arg3__)
                                    {
#line hidden
                                        Release__.Apply(flagQubit);
                                    }
                                }
                            }
                        }
#line hidden
                        catch
                        {
                            __arg2__ = false;
                            throw;
                        }
#line hidden
                        finally
                        {
                            if (__arg2__)
                            {
#line hidden
                                Release__.Apply(referenceRegister);
                            }
                        }
                    }
                }
#line hidden
                catch
                {
                    __arg1__ = false;
                    throw;
                }
#line hidden
                finally
                {
                    if (__arg1__)
                    {
#line hidden
                        Release__.Apply(dataRegister);
                    }
                }
            }

#line hidden
            return QVoid.Instance;
        }

        ;
        public override void __Init__()
        {
            this.Allocate__ = this.__Factory__.Get<Allocate>(typeof(global::Microsoft.Quantum.Intrinsic.Allocate));
            this.Release__ = this.__Factory__.Get<Release>(typeof(global::Microsoft.Quantum.Intrinsic.Release));
            this.Microsoft__Quantum__Arrays__Zipped = this.__Factory__.Get<ICallable>(typeof(global::Microsoft.Quantum.Arrays.Zipped<,>));
            this.Microsoft__Quantum__Intrinsic__X = this.__Factory__.Get<IUnitary<Qubit>>(typeof(global::Microsoft.Quantum.Intrinsic.X));
            this.AmplitudeTransduction__UnifPrime = this.__Factory__.Get<IUnitary<(IQArray<Qubit>,IQArray<Qubit>)>>(typeof(global::AmplitudeTransduction.UnifPrime));
            this.ResetAll__ = this.__Factory__.Get<ICallable<IQArray<Qubit>, QVoid>>(typeof(global::Microsoft.Quantum.Intrinsic.ResetAll));
            this.Reset__ = this.__Factory__.Get<ICallable<Qubit, QVoid>>(typeof(global::Microsoft.Quantum.Intrinsic.Reset));
        }

        public override IApplyData __DataIn__(QVoid data) => data;
        public override IApplyData __DataOut__(QVoid data) => data;
        public static System.Threading.Tasks.Task<QVoid> Run(IOperationFactory __m__)
        {
            return __m__.Run<TestUnifPrime, QVoid, QVoid>(QVoid.Instance);
        }
    }
}