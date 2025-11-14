using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEditor;
using UnityEngine;
using UnityEngine.UIElements;
public class PhysicsHandler : MonoBehaviour
{
    [Header("Inertia")]
    [SerializeField] private float Ixx = 100;
    [SerializeField] private float Iyy = 100;
    [SerializeField] private float Izz = 100;

    [SerializeField] private float Ixy = 0;
    [SerializeField] private float Ixz = 0;
    [SerializeField] private float Iyz = 0;

    [Header("Rotations")]
    [SerializeField] public Vector3 RPY;
    [SerializeField] private Vector3 bodyRot;

    [Header("Mass")]
    [SerializeField] public float mass = 600;
    [SerializeField] public float g = 9.81f;

    [Header("Linear")]
    [SerializeField] public Vector3 XYZ;
    [SerializeField] public Vector3 bodyVel;

    [HideInInspector] public Quaternion quatRH;
    private Quaternion quatLH;
    private Vector3 pqr;           // Current Body Rotation for 3D Rigidbody Rotation Calculations (RADS)
    private Vector3 pqrpast;       // Past Body Rotation Rates for 3D Rigidbody Rotation Calculations
    private Vector3 pqrdot;
    private Vector3 pqrdotpast;
    private Vector3 RPYint;

    private Vector3 uvw;
    private Vector3 uvwdot;
    private Vector3 XYZdot;
    private Vector3 visXYZ;

    [HideInInspector] public Vector3 F;
    [HideInInspector] public Vector3 M;
    [HideInInspector] public Vector3 P;
    [HideInInspector] public Vector3 H;

    private float DebugLogTimer = 0f;

    void Start()
    {
        start_toQuaternion();
        start_variableDef();
        transform.SetLocalPositionAndRotation(visXYZ, quatLH);
    }


    void FixedUpdate()
    {
        CalculateAccel();
        CalculateVel();
        CalculatePos();
        UpdateFlightData();
    }

    void CalculateAccel()
    {
        CalcAngAcc();
        CalcLatAcc();
    }

    void CalcAngAcc()
    {
        pqrdot.x = ((Iyy - Izz) * pqrpast.y * pqrpast.z + Ixz * (pqrdotpast.z + pqrpast.x * pqrpast.y) + M.x) / Ixx;
        pqrdot.y = ((Izz - Ixx) * pqrpast.x * pqrpast.z + Ixz * (pqrpast.z * pqrpast.z - pqrpast.x * pqrpast.x) + M.y) / Iyy;
        pqrdot.z = ((Ixx - Iyy) * pqrpast.x * pqrpast.y + Ixz * (pqrdotpast.x + pqrpast.y * pqrpast.z) + M.z) / Izz;

        pqrdotpast = pqrdot;
    }

    void CalcLatAcc()
    {
        (float sp, float cp, float sr, float cr) = EulerDet();
        uvwdot.x = (F.x / mass - g * sp) + pqrpast.z * uvw.y - pqrpast.y * uvw.z;
        uvwdot.y = (F.y / mass + g * cp * sr) + pqrpast.x * uvw.z - pqrpast.z * uvw.x;
        uvwdot.z = (F.z / mass + g * cp * cr) + pqrpast.y * uvw.x - pqrpast.x * uvw.y;

    }

    (float, float, float, float) EulerDet()
    {
        float sp = Mathf.Sin(RPY.y * Mathf.Deg2Rad);
        float cp = Mathf.Cos(RPY.y * Mathf.Deg2Rad);
        float sr = Mathf.Sin(RPY.x * Mathf.Deg2Rad);
        float cr = Mathf.Cos(RPY.x * Mathf.Deg2Rad);

        return (sp, cp, sr, cr);
    }
    void CalculateVel()
    {

        CalcAngVel();
        CalcLatVel();

    }

    void CalcAngVel()
    {

        pqrpast = pqr;
        pqr += pqrdot * Time.fixedDeltaTime;
        H = new Vector3(Ixx * pqr.x, Iyy * pqr.y, Izz * pqr.z);
    }

    void CalcLatVel()
    {
        uvw += uvwdot * Time.fixedDeltaTime;

        float cr = Mathf.Cos(RPYint.x);
        float cp = Mathf.Cos(RPYint.y);
        float cy = Mathf.Cos(RPYint.z);

        float sr = Mathf.Sin(RPYint.x);
        float sp = Mathf.Sin(RPYint.y);
        float sy = Mathf.Sin(RPYint.z);

        XYZdot.x = cp * cy * uvw.x + (sr * sp * cy - cr * sy) * uvw.y + (cr * sp * cy + sr * sy) * uvw.z;
        XYZdot.y = cp * sy * uvw.x + (sr * sp * sy + cr * cy) * uvw.y + (cr * sp * sy - sr * cy) * uvw.z;
        XYZdot.z = (-sp) * uvw.x + (sr * cp) * uvw.y + (cr * cp) * uvw.z;
    }
    void CalculatePos()
    {
        CalcAngPos();
        CalcLatPos();
    }

    void CalcAngPos()
    {
        Matrix4x4 Omega = new Matrix4x4();
        Omega = CreateOmega(Omega);
        Vector4 qInter = new Vector4(quatRH.x, quatRH.y, quatRH.z, quatRH.w);
        qInter = Omega * qInter;

        quatRH = new Quaternion(qInter.x, qInter.y, qInter.z, qInter.w);
        quatRH = quatRH.normalized;
        quatLH = new Quaternion(quatRH.z, quatRH.x, quatRH.y, -quatRH.w);
    }

    void CalcLatPos()
    {
        XYZ += XYZdot * Time.fixedDeltaTime;
        visXYZ = new Vector3(XYZ.x, -XYZ.z, -XYZ.y);
    }
    Matrix4x4 CreateOmega(Matrix4x4 omega)
    {
        Vector3 pqrM = pqr;

        omega.SetRow(0, (new Vector4(0, -pqrM.x, -pqrM.y, -pqrM.z)) * 0.5f * Time.fixedDeltaTime);
        omega.SetRow(1, (new Vector4(pqrM.x, 0, pqrM.z, -pqrM.y)) * 0.5f * Time.fixedDeltaTime);
        omega.SetRow(2, (new Vector4(pqrM.y, -pqrM.z, 0, pqrM.x)) * 0.5f * Time.fixedDeltaTime);
        omega.SetRow(3, (new Vector4(pqrM.z, pqrM.y, -pqrM.x, 0)) * 0.5f * Time.fixedDeltaTime);

        omega = matAdd(omega, Matrix4x4.identity);

        return omega;
    }

    Matrix4x4 matAdd(Matrix4x4 mat1, Matrix4x4 mat2)
    {
        Matrix4x4 mat3 = new Matrix4x4();
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                mat3[i, j] = mat1[i, j] + mat2[i, j];
            }
        }

        return mat3;
    }

    void UpdateFlightData()
    {
        transform.SetLocalPositionAndRotation(visXYZ, quatLH);
        bodyRot = pqr * Mathf.Rad2Deg;
        bodyVel = uvw;
        toEuler();
        RPY = RPYint * Mathf.Rad2Deg;
    }
    public void start_toQuaternion()
    {
        float hcr = Mathf.Cos(Mathf.Deg2Rad * 0.5f * -RPY.z);
        float hsr = Mathf.Sin(Mathf.Deg2Rad * 0.5f * -RPY.z);
        float hcp = Mathf.Cos(Mathf.Deg2Rad * 0.5f * -RPY.y);
        float hsp = Mathf.Sin(Mathf.Deg2Rad * 0.5f * -RPY.y);
        float hcy = Mathf.Cos(Mathf.Deg2Rad * 0.5f * RPY.x);
        float hsy = Mathf.Sin(Mathf.Deg2Rad * 0.5f * RPY.x);

        quatRH.w = hcr * hcp * hcy + hsr * hsp * hsy;
        quatRH.x = hsr * hcp * hcy - hcr * hsp * hsy;
        quatRH.y = hcr * hsp * hcy + hsr * hcp * hsy;
        quatRH.z = hcr * hcp * hsy - hsr * hsp * hcy;

        quatLH = new Quaternion(-quatRH.x, quatRH.z, quatRH.y, quatRH.w);
    }

    public void start_variableDef()
    {
        pqr = bodyRot * Mathf.Deg2Rad;
        pqrpast = pqr;
        uvw = bodyVel;
        RPYint = RPY * Mathf.Deg2Rad;
        visXYZ = new Vector3(XYZ.x, -XYZ.z, -XYZ.y);
    }

    void toEuler()
    {
        float srcp = 2 * (quatRH.w * quatRH.x + quatRH.y * quatRH.z);
        float crcp = 1 - 2 * (quatRH.x * quatRH.x + quatRH.y * quatRH.y);
        RPYint.z = -Mathf.Atan2(srcp, crcp);

        float sp = Mathf.Sqrt(1 + 2 * (quatRH.w * quatRH.y - quatRH.x * quatRH.z));
        float cp = Mathf.Sqrt(1 - 2 * (quatRH.w * quatRH.y - quatRH.x * quatRH.z));
        RPYint.y = -(2 * Mathf.Atan2(sp, cp) - Mathf.PI / 2);
        if (float.IsNaN(RPYint.y))
        {
            RPYint.y = 0;
        }

        float sycp = 2 * (quatRH.w * quatRH.z + quatRH.x * quatRH.y);
        float cycp = 1 - 2 * (quatRH.y * quatRH.y + quatRH.z * quatRH.z);
        RPYint.x = Mathf.Atan2(sycp, cycp);
    }
}
