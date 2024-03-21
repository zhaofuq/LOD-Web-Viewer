import localforage from "localforage";


export async function canInsetToLocalDB() {
	let flag = true;
	if (navigator.storage && navigator.storage.estimate) {
		const quota = await navigator.storage.estimate();
		const percentageUsed = quota.usage / quota.quota;
		console.log('DB使用:', percentageUsed)
		if (percentageUsed > 0.8) {
			localforage.clear()
			flag = false;
		}
	}
	return flag
}

export async function localforageSetItem(key, value) {
	let canInsert = canInsetToLocalDB()
	if (!canInsert) {
		await localforage.clear()
		localforage.setItem(key, value);
	}
	else {
		localforage.setItem(key, value);
	}
}
export async function localforageGetItem(key) {
	const localData = await localforage.getItem(key);
	return localData;
}