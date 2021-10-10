extern "C" DLLEXPORT bool SKSEAPI SKSEPlugin_Query(const SKSE::QueryInterface* a_skse, SKSE::PluginInfo* a_info)
{
#ifndef NDEBUG
	auto sink = std::make_shared<spdlog::sinks::msvc_sink_mt>();
#else
	auto path = logger::log_directory();
	if (!path) {
		return false;
	}

	*path /= Version::PROJECT;
	*path += ".log"sv;
	auto sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(path->string(), true);
#endif

	auto log = std::make_shared<spdlog::logger>("global log"s, std::move(sink));

#ifndef NDEBUG
	log->set_level(spdlog::level::trace);
#else
	log->set_level(spdlog::level::info);
	log->flush_on(spdlog::level::warn);
#endif

	spdlog::set_default_logger(std::move(log));
	spdlog::set_pattern("%g(%#): [%^%l%$] %v"s);

	logger::info(FMT_STRING("{} v{}"), Version::PROJECT, Version::NAME);

	a_info->infoVersion = SKSE::PluginInfo::kVersion;
	a_info->name = Version::PROJECT.data();
	a_info->version = Version::MAJOR;

	if (a_skse->IsEditor()) {
		logger::critical("Loaded in editor, marking as incompatible"sv);
		return false;
	}

	const auto ver = a_skse->RuntimeVersion();
	if (ver < SKSE::RUNTIME_1_5_39) {
		logger::critical(FMT_STRING("Unsupported runtime version {}"), ver.string());
		return false;
	}

	return true;
}

namespace StaminaCosts
{
	const float STAMINA_PER_LVL = 10.0f;
	const size_t A_shield_size = 9;
	const float A_shield1[A_shield_size] = {
		0.01785038659f, -0.00003066918536f, -0.7363035728f,
		-0.0106493435f, -0.0008806043337f, 0.01624281144f,
		0.3140566047f, -0.471075562f, 10.84342129f
	};
	const float A_shield2[A_shield_size] = {
		0.01178125515f, -0.00002024166234f, -0.485960358f,
		-0.007028566711f, -0.0005811988603f, 0.01072025555f,
		0.2072773591f, -0.3109098709f, 7.156658054f
	};
	const size_t A_sz = 3;

	const float A1[] = {
		-0.0001200869909f, 0.00280348075f, -0.01826128321f,
		0.02111060035f, -0.4183658809f, 1.970622115f,
		-0.02099051336f, 0.4155624002f, -1.952360832f,
		0.000075112713f, -0.001709687183f, 0.01066609571f,
		-0.01282553565f, 0.231304605f, -0.5246021088f,
		0.01275042293f, -0.2295949179f, 0.5139360131f,
		-0.00001095903224f, 0.0002442156723f, -0.001416952899f,
		0.00186941694f, -0.0324241964f, 0.03899699052f,
		-0.002958330682f, -0.03906864608f, 5.034768462f
	};
	const float A2[] = {
		-0.00005028357888f, 0.001152142022f, -0.005540211826f,
		0.009733653994f, -0.1968288162f, 0.4586011771f,
		-0.009683370415f, 0.1956766742f, -0.4530609653f,
		0.00003238059428f, -0.0007098323266f, 0.002767430244f,
		-0.005793170465f, 0.1063998537f, 0.04773309646f,
		0.00576078987f, -0.1056900214f, -0.0505005267f,
		-0.000004728615976f, 0.00009946901221f, -0.000284950707f,
		0.000795935801f, -0.0137692291f, -0.03263028013f,
		-0.000791207185f, 0.01366976009f, 1.032915231f
	};

	const float A_jump[] = { 0.4308323564f, 13.28057835f };
	const float A_regn[] = { 0.6462485346f, 9.920867526f };
	const float A_bwdr[] = { -1.0f, 25.0f };
	const float A_Kmul[] = { 6.787330317f, 0.2488687783f, 0.1787330317f };

	template <size_t N>
	float mulN(const float* a, const float* b)
	{
		float ans = 0.0f;
		for (size_t i = 0; i < N; i++) {
			ans += a[i] * b[i];
		}
		return ans;
	}

	float mul(const float* a, const float* b)
	{
		return mulN<A_sz>(a, b);
	}

	float mul2(const float* a, const float* b)
	{
		return mulN<2>(a, b);
	}

	void mull_VM_T(const float* v, const float* M_T, float* ans)
	{
		for (size_t i = 0; i < A_sz; i++) {
			ans[i] = mul(v, &M_T[i * A_sz]);
		}
	}

	void mull_vA3(const float* w, const float* A, float* wA)
	{
		for (size_t i = 0; i < A_sz; i++) {
			mull_VM_T(w, &A[i * A_sz * A_sz], &wA[i * A_sz]);
		}
	}

	void mull_Av(const float* A, const float* w, float* Aw)
	{
		for (size_t i = 0; i < A_sz; i++) {
			Aw[i] = mul(&A[i * A_sz], w);
		}
	}

	constexpr float MIN_K = 0.3148688047f;
	constexpr float MID_K = 0.6798928223f;
	constexpr float MAX_K = 0.8840841683f;

	float kernel(float _k)
	{
		constexpr float Ax = 0.35f, k = 11.66180758f;
		if (_k <= Ax)
			return _k * _k * _k * k;
		if (_k > 1.0f)
			_k = 0.99f;
		constexpr float m = 5.512008447f, a = 5.571428571f;
		return 1.0f - m * std::pow(1 - _k, a);
	}

	float get_k(float l, float ST)
	{
		if (l > 1)
			return kernel((ST - 100.0f) / ((l - 1.0f) * STAMINA_PER_LVL));
		else
			return MID_K;
	}

	float getAttackCount(float l, float ST, float w, bool ispower, bool isPlayer)
	{
		float k = isPlayer ? get_k(l, ST) : MID_K;

		float W[3] = { w * w, w, 1 };
		float L[3] = { l * l, l, 1 };
		float K[3] = { k * k, k, 1 };
		const float* M = A1;
		if (ispower)
			M = A2;

		float res1[9];
		float res2[3];
		mull_vA3(W, M, res1);
		mull_Av(res1, L, res2);
		return mul(res2, K);
	}

	float getShieldCost_(float l, float ST, float w_weap, float w_shield, bool ispower)
	{
		const float data[A_shield_size] = { w_weap * w_weap, l * l, w_weap, l, l * w_weap, ST, log(l), w_shield, 1 };
		const float* A = A_shield1;
		if (ispower)
			A = A_shield2;
		auto F = mulN<9>(data, A);
		return F < 0.000001f ? 5000.0f : ST / F;
	}

	float safeGetWeaponWeight(RE::Actor* a, bool left)
	{
		float w = 1.0;
		auto weap = a->GetEquippedObject(left);
		if (weap)
			w = weap->GetWeight();
		return w;
	}

	float getShieldCost(RE::StaticFunctionTag*, RE::Actor* target, RE::Actor* heitor, bool ispower, bool isPlayer)
	{
		float l = target->GetLevel();
		float ST = isPlayer ? target->GetBaseActorValue(RE::ActorValue::kStamina) : (l - 1) * 10 / 3.0f + 100;
		return getShieldCost_(l, ST, safeGetWeaponWeight(heitor, false), safeGetWeaponWeight(target, true), ispower) / 2.0f;
	}

	float safeGetArmorWeight(RE::Actor* a, uint32_t slot)
	{
		auto changes = static_cast<RE::ExtraContainerChanges*>(a->extraList.GetByType(RE::ExtraDataType::kContainerChanges));
		if (changes) {
			auto armo = changes->changes->GetArmorInSlot(slot + 30);
			if (armo) {
				return armo->GetWeight();
			}
		}
		return 0.0f;
	}

	float getCuriassWeight(RE::Actor* a)
	{
		return safeGetArmorWeight(a, 0x2);
	}

	float getArmopartsWeight(RE::Actor* a)
	{
		return safeGetArmorWeight(a, 0x1) + safeGetArmorWeight(a, 0x3) + safeGetArmorWeight(a, 0x7);
	}

	float get_ArmorModifier(RE::Actor* a)
	{
		float x = getCuriassWeight(a), y = getArmopartsWeight(a);
		const float alpha = 0.2f;
		const float k = 0.3f, b = 0.33f;
		return b + k * (log(x + 1.0f) + alpha * log(y + 1.0f)) / (1.0f + alpha);
	}

	float getAttackCost(RE::StaticFunctionTag*, RE::Actor* a, bool left, bool isPlayer)
	{
		float w = safeGetWeaponWeight(a, left), l = a->GetLevel(),
			  ST = a->GetBaseActorValue(RE::ActorValue::kStamina),
			  F1 = getAttackCount(l, ST, w, false, isPlayer);
		return ST * get_ArmorModifier(a) / F1;
	}

	float getAttackDelta(RE::StaticFunctionTag*, RE::Actor* a, bool left, bool isPlayer)
	{
		float w = safeGetWeaponWeight(a, left), l = a->GetLevel(),
			  ST = a->GetBaseActorValue(RE::ActorValue::kStamina),
			  F1 = getAttackCount(l, ST, w, false, isPlayer), F2 = getAttackCount(l, ST, w, true, isPlayer);
		return (ST / F2 - ST * 3 / F1) * get_ArmorModifier(a);
	}

	float safeGetWeight(RE::Actor* a, RE::BGSBipedObjectForm::BipedObjectSlot slot)
	{
		auto armor = a->GetWornArmor(slot);
		if (!armor)
			return 0.0f;
		else
			return armor->GetWeight();
	}

	float safeGetWeight(RE::Actor* a, bool left = false)
	{
		auto weapon = a->GetEquippedObject(left);
		if (!weapon)
			return 0.0f;
		else
			return weapon->GetWeight();
	}

	float getCarryWeight(RE::Actor* a, bool withweapons)
	{
		// bug -- dont count diadems
		float ans = 0.0f;
		ans += safeGetWeight(a, RE::BGSBipedObjectForm::BipedObjectSlot::kHead);
		ans += safeGetWeight(a, RE::BGSBipedObjectForm::BipedObjectSlot::kBody);
		ans += safeGetWeight(a, RE::BGSBipedObjectForm::BipedObjectSlot::kHands);
		ans += safeGetWeight(a, RE::BGSBipedObjectForm::BipedObjectSlot::kFeet);
		if (withweapons) {
			// bug -- dont count left hand
			//ans += safeGetWeight(a, true);
			ans += safeGetWeight(a, false);
		}
		return ans;
	}

	float getJumpCost(RE::StaticFunctionTag*, RE::Actor* a)
	{
		auto cw = getCarryWeight(a, true);
		const float b[] = { cw, 1.0f };
		return mul2(b, A_jump);
	}

	float getNewStaminaRate(RE::StaticFunctionTag*, RE::Actor* a)
	{
		auto cw = getCarryWeight(a, true);
		const float b[] = { cw, 1.0f };
		auto x = mul2(b, A_regn);
		auto pers = a->GetActorValue(RE::ActorValue::kHealth) / a->GetBaseActorValue(RE::ActorValue::kHealth);
		return 100.0f * pers * pers / x;
	}

	float getBowCost(RE::StaticFunctionTag*, RE::Actor* a, bool isPlayer)
	{
		auto bow = a->GetEquippedObject(false);
		if (!bow)
			return 0.0f;
		else {
			float k = isPlayer ? get_k(a->GetLevel(), a->GetBaseActorValue(RE::ActorValue::kStamina)) : MID_K;

			auto w = bow->GetWeight();
			const float b[] = { w, 1.0f };
			float attacks = mul2(b, A_bwdr);
			float K[] = { k * k, k, 1 };
			attacks = std::max(1.0f, attacks) * mul(K, A_Kmul);
			return 200.0f / attacks;
		}
	}
}

namespace Staggers
{
	// between t's back and c's attack
	float getAngle(float t, float c)
	{
		auto ans = 180.0f * abs(c - t) / 3.1415926f;
		if (ans > 180.0f)
			ans = 360.0f - ans;
		return ans;
	}

	float scale(float x, float oldL, float oldR, float newL, float newR)
	{
		return (x - oldL) * (newR - newL) / (oldR - oldL) + newL;
	}

	bool iscreature(RE::Actor* a)
	{
		return !a->GetRace() || (a->GetRace()->validEquipTypes.underlying() & 1) == 0;
	}

	float getEquippementInfo_armo(RE::Actor* target)
	{
		if (iscreature(target)) {
			return target->GetBaseActorValue(RE::ActorValue::kHealth) * (21.0f / 200.0f);
		} else {
			auto armo = target->GetWornArmor(RE::BGSBipedObjectForm::BipedObjectSlot::kBody);
			return armo ? armo->GetWeight() : 0.0f;
		}
	}

	float getEquippementInfo_weap(RE::Actor* heitor)
	{
		if (iscreature(heitor)) {
			return heitor->GetBaseActorValue(RE::ActorValue::kStamina) * (18.0f / 200.0f);
		} else {
			auto weap = heitor->GetEquippedObject(false);
			if (!weap)
				weap = heitor->GetEquippedObject(true);
			return weap ? weap->GetWeight() : 0.0f;
		}
	}

	std::pair<float, float> getEquippementInfo(RE::Actor* target, RE::Actor* heitor)
	{
		return { getEquippementInfo_armo(target), getEquippementInfo_weap(heitor) };
	}

	void setParams(RE::Actor* target, RE::Actor* heitor, bool ispower, float* params)
	{
		auto info = getEquippementInfo(target, heitor);
		params[0] = std::min(360.0f, scale(info.first, 8.0f, 34.0f, 180.0f, 90.0f));
		params[2] = std::max(0.0f, scale(info.second, 14.0f, 25.0f, 10.0f, 15.0f));
		params[1] = params[2] * 0.3f;

		if (ispower) {
			params[0] = 360.0f;
			params[1] = params[2] * 0.5f;
		}

		float delta = heitor->GetBaseActorValue(RE::ActorValue::kHealth) - target->GetBaseActorValue(RE::ActorValue::kHealth);
		if (delta > 0.0f) {
			float needdelta = (360.0f - params[0]) * (100.0f / 90.0f);
			if (needdelta < delta) {
				delta -= needdelta;
				params[0] = 360.0f;
				needdelta = (params[2] - params[1]) * (25.0f / 1.0f);
				if (needdelta < delta) {
					params[1] = params[2];
					delta -= needdelta;
					params[2] += (1.0f / 25.0f) * delta;
				} else {
					params[1] += (1.0f / 25.0f) * delta;
				}
			} else {
				params[0] += (90.0f / 100.0f) * delta;
			}
		}
		params[1] *= 0.2f, params[2] *= 0.2f;
		params[2] = std::min(4.0f, params[2]);
	}

	float getStaggerTime_(RE::Actor* target, RE::Actor* heitor, bool ispower)
	{
		float angle = getAngle(target->GetAngleZ(), heitor->GetAngleZ());
		float params[3];
		setParams(target, heitor, ispower, params);
		float maxAngle = params[0] / 2.0f;
		if (angle > maxAngle || params[2] <= 0.0f)
			return -1.0f;
		return std::min(4.0f, scale(angle, maxAngle, 0.0f, params[1], params[2]));
	}

	float getStaggerTime(RE::StaticFunctionTag*, RE::Actor* target, RE::Actor* heitor, bool ispower)
	{
		return getStaggerTime_(target, heitor, ispower);
	}

	int getStaggerType(RE::StaticFunctionTag*, RE::Actor* target, RE::Actor* heitor, bool ispower, float mult)
	{
		float time = getStaggerTime_(target, heitor, ispower) * mult;
		if (time >= 3.0f)
			return 2;
		if (time >= 1.0f)
			return 1;
		if (time >= 0)
			return 0;
		return -1;
	}
}

float FenixLog(RE::StaticFunctionTag*, float x)
{
	return log(x);
}

bool RegisterFuncs(RE::BSScript::IVirtualMachine* a_vm)
{
	a_vm->RegisterFunction("FenixLog", "f314IM_Utility", FenixLog);

	a_vm->RegisterFunction("getAttackCost", "f314IM_Utility", StaminaCosts::getAttackCost);
	a_vm->RegisterFunction("getAttackDelta", "f314IM_Utility", StaminaCosts::getAttackDelta);
	a_vm->RegisterFunction("getShieldCost", "f314IM_Utility", StaminaCosts::getShieldCost);
	a_vm->RegisterFunction("getJumpCost", "f314IM_Utility", StaminaCosts::getJumpCost);
	a_vm->RegisterFunction("getNewStaminaRate", "f314IM_Utility", StaminaCosts::getNewStaminaRate);
	a_vm->RegisterFunction("getBowCost", "f314IM_Utility", StaminaCosts::getBowCost);

	a_vm->RegisterFunction("getStaggerType", "f314IM_Utility", Staggers::getStaggerType);
	a_vm->RegisterFunction("getStaggerTime", "f314IM_Utility", Staggers::getStaggerTime);
	return true;
}

extern "C" DLLEXPORT bool SKSEAPI SKSEPlugin_Load(const SKSE::LoadInterface* a_skse)
{
	SKSE::Init(a_skse);

	auto papyrus = SKSE::GetPapyrusInterface();
	if (!papyrus->Register(RegisterFuncs)) {
		return false;
	}

	logger::info("loaded");

	return true;
}
