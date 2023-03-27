import json
import KratosMultiphysics

#j = json.loads('{"a" : 1,"b" : { "c" : 2},"d" : 3}')

a = KratosMultiphysics.Parameters("""
    {
        "a" : 1, //pippo
        "b" : { "c" : 2},
        "d" : 3
    }
    """)

print(a)
print(a["a"].GetInt())
print(a[jsonPointer("/b/c")].GetInt())