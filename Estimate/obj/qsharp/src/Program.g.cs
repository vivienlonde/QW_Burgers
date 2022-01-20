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

[assembly: CallableDeclaration("{\"Kind\":{\"Case\":\"Operation\"},\"QualifiedName\":{\"Namespace\":\"Estimate\",\"Name\":\"TestAmplitudeTransduction\"},\"Attributes\":[{\"TypeId\":{\"Case\":\"Value\",\"Fields\":[{\"Namespace\":\"Microsoft.Quantum.Core\",\"Name\":\"EntryPoint\",\"Range\":{\"Case\":\"Value\",\"Fields\":[{\"Item1\":{\"Line\":1,\"Column\":2},\"Item2\":{\"Line\":1,\"Column\":12}}]}}]},\"TypeIdRange\":{\"Case\":\"Value\",\"Fields\":[{\"Item1\":{\"Line\":1,\"Column\":2},\"Item2\":{\"Line\":1,\"Column\":12}}]},\"Argument\":{\"Item1\":{\"Case\":\"UnitValue\"},\"Item2\":[],\"Item3\":{\"Case\":\"UnitType\"},\"Item4\":{\"IsMutable\":false,\"HasLocalQuantumDependency\":false},\"Item5\":{\"Case\":\"Value\",\"Fields\":[{\"Item1\":{\"Line\":1,\"Column\":12},\"Item2\":{\"Line\":1,\"Column\":14}}]}},\"Offset\":{\"Item1\":16,\"Item2\":4},\"Comments\":{\"OpeningComments\":[],\"ClosingComments\":[]}}],\"Modifiers\":{\"Access\":{\"Case\":\"DefaultAccess\"}},\"SourceFile\":\"C:\\\\Users\\\\vivie\\\\OneDrive\\\\Documents\\\\GitHub\\\\QW_Burgers\\\\Estimate\\\\Program.qs\",\"Position\":{\"Item1\":17,\"Item2\":4},\"SymbolRange\":{\"Item1\":{\"Line\":1,\"Column\":11},\"Item2\":{\"Line\":1,\"Column\":36}},\"ArgumentTuple\":{\"Case\":\"QsTuple\",\"Fields\":[[]]},\"Signature\":{\"TypeParameters\":[],\"ArgumentType\":{\"Case\":\"UnitType\"},\"ReturnType\":{\"Case\":\"UnitType\"},\"Information\":{\"Characteristics\":{\"Case\":\"EmptySet\"},\"InferredInformation\":{\"IsSelfAdjoint\":false,\"IsIntrinsic\":false}}},\"Documentation\":[]}")]
[assembly: SpecializationDeclaration("{\"Kind\":{\"Case\":\"QsBody\"},\"TypeArguments\":{\"Case\":\"Null\"},\"Information\":{\"Characteristics\":{\"Case\":\"EmptySet\"},\"InferredInformation\":{\"IsSelfAdjoint\":false,\"IsIntrinsic\":false}},\"Parent\":{\"Namespace\":\"Estimate\",\"Name\":\"TestAmplitudeTransduction\"},\"Attributes\":[],\"SourceFile\":\"C:\\\\Users\\\\vivie\\\\OneDrive\\\\Documents\\\\GitHub\\\\QW_Burgers\\\\Estimate\\\\Program.qs\",\"Position\":{\"Item1\":17,\"Item2\":4},\"HeaderRange\":{\"Item1\":{\"Line\":1,\"Column\":11},\"Item2\":{\"Line\":1,\"Column\":36}},\"Documentation\":[]}")]
#line hidden
namespace Estimate
{
    [SourceLocation("C:\\Users\\vivie\\OneDrive\\Documents\\GitHub\\QW_Burgers\\Estimate\\Program.qs", OperationFunctor.Body, 18, -1)]
    public partial class TestAmplitudeTransduction : Operation<QVoid, QVoid>, ICallable
    {
        public TestAmplitudeTransduction(IOperationFactory m) : base(m)
        {
        }

        String ICallable.Name => "TestAmplitudeTransduction";
        String ICallable.FullName => "Estimate.TestAmplitudeTransduction";
        public static EntryPointInfo<QVoid, QVoid> Info => new EntryPointInfo<QVoid, QVoid>(typeof(TestAmplitudeTransduction));
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

        protected ICallable<(Int64,IQArray<Qubit>), Microsoft.Quantum.Arithmetic.FixedPoint> Microsoft__Quantum__Arithmetic__FixedPoint
        {
            get;
            set;
        }

        protected IUnitary<Microsoft.Quantum.Arithmetic.FixedPoint> AmplitudeTransduction__AmplitudeTransduction
        {
            get;
            set;
        }

        public override Func<QVoid, QVoid> __Body__ => (__in__) =>
        {
#line 19 "C:\\Users\\vivie\\OneDrive\\Documents\\GitHub\\QW_Burgers\\Estimate\\Program.qs"
            var nbQubits = 5L;
#line hidden
            {
#line 20 "C:\\Users\\vivie\\OneDrive\\Documents\\GitHub\\QW_Burgers\\Estimate\\Program.qs"
                var transitionProbabilityQubitRegister = Allocate__.Apply(nbQubits);
#line hidden
                bool __arg1__ = true;
                try
                {
#line 21 "C:\\Users\\vivie\\OneDrive\\Documents\\GitHub\\QW_Burgers\\Estimate\\Program.qs"
                    var transitionProbability = new Microsoft.Quantum.Arithmetic.FixedPoint((0L, transitionProbabilityQubitRegister));
#line 22 "C:\\Users\\vivie\\OneDrive\\Documents\\GitHub\\QW_Burgers\\Estimate\\Program.qs"
                    AmplitudeTransduction__AmplitudeTransduction.Apply(transitionProbability);
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
                        Release__.Apply(transitionProbabilityQubitRegister);
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
            this.Microsoft__Quantum__Arithmetic__FixedPoint = this.__Factory__.Get<ICallable<(Int64,IQArray<Qubit>), Microsoft.Quantum.Arithmetic.FixedPoint>>(typeof(global::Microsoft.Quantum.Arithmetic.FixedPoint));
            this.AmplitudeTransduction__AmplitudeTransduction = this.__Factory__.Get<IUnitary<Microsoft.Quantum.Arithmetic.FixedPoint>>(typeof(global::AmplitudeTransduction.AmplitudeTransduction));
        }

        public override IApplyData __DataIn__(QVoid data) => data;
        public override IApplyData __DataOut__(QVoid data) => data;
        public static System.Threading.Tasks.Task<QVoid> Run(IOperationFactory __m__)
        {
            return __m__.Run<TestAmplitudeTransduction, QVoid, QVoid>(QVoid.Instance);
        }
    }
}