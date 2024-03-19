
import { WorkerPool } from "./gaussian-utils/WorkerPool.js";
import { OctreeLoader } from "./gaussian-utils/gaussian-octree-loader.js";

export const workerPool = new WorkerPool();

export const version = {
	major: 1,
	minor: 8,
	suffix: '.0'
};

// export let pointBudget = 1 * 1000 * 1000;
// export let framenumber = 0;
// export let numNodesLoading = 0;
// export let maxNodesLoading = 4;

export const potree = {
	pointBudget: 1 * 1000 * 1000,
	framenumber: 0,
	numNodesLoading: 0,
	maxNodesLoading: 4
}

export async function loadPointCloud(path, name, octreeFileUrl, callback) {

	let loaded = function (e) {
		e.geometry.name = name;
		callback(e);
	};

	let promise = new Promise(async (resolve, reject) => {
		// load pointcloud
		if (!path) {
			reject(new Error("Path is undefined or empty"));
		} else if (path.indexOf('metadata.json') > 0) {
			OctreeLoader.load(path, octreeFileUrl).then(e => {
				let geometry = e.geometry;
				if (!geometry) {
					console.error(new Error(`failed to load point cloud from URL: ${path}`));
				} else {

					resolve({ type: 'geometry', geometry: geometry });
				}
			}).catch(e => {
				console.error(new Error(`failed to load point cloud from URL: ${path}`));
				reject(e);
			});
		} else {
			//callback({'type': 'loading_failed'});
			reject(new Error(`Invalid path for point cloud: ${path}`));
		}
	});

	if (callback) {
		promise.then(pointcloud => {
			loaded(pointcloud);
		});
	} else {
		return promise;
	}
}